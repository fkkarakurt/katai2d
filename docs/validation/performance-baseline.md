# Performans kampanyası — taban çizgisi ve adım adım kazanımlar

Amaç: hesap motorunu doğruluğu **birebir koruyarak** hızlandırmak. Yöntem: önce ölç
(`study_perf`, GUI compute yolu üzerinden gerçek senaryolar), sonra en büyük kalemden
başlayarak tek tek optimize et; her adımda (a) `study_perf` hızlanmayı, (b) tam ctest
doğruluğun korunduğunu göstermek zorunda.

Çalıştırma:

```
cmake --build build/msvc-rwdi --target study_perf
build/msvc-rwdi/bin/study_perf.exe
```

Enstrümantasyon: `NewtonResult::Timings` (nonlinear_solver içi, her zaman açık;
iterasyon granülaritesinde chrono → maliyeti ihmal edilebilir). GUI'ye
`SolveResult.timings/iterations` ile taşınır (PLAXIS'in calculation-time raporu gibi).

Makine: Win11, MSVC RelWithDebInfo, MKL 2026.0 sequential statik (tek iş parçacığı).

## Senaryolar

| Senaryo | Ne ölçer |
|---|---|
| LE 40×20 gravity (coarse/fine) | Lineer çözüm baskın yol (1 adım, ~2 iterasyon) |
| MC footing collapse (Reissner φ=20, 1.8×q_ult) | Çok iterasyonlu artımsal limit analizi + cutback |
| HS footing service (150 kPa, K0+SumMstage) | Bünye entegrasyonu (substepping) ağır yol, 40 adım |

## Taban çizgisi (commit 3cc15a8 + enstrümantasyon, 2026-06-12)

```
scenario                    elems    dofs iters | tan-asm res-asm    csr  pardiso  other |    total
LE 40x20 gravity coarse      2333    9576     2 |    0.03    0.00   0.01     0.14   0.00 |     0.18 s
LE 40x20 gravity fine       11926   48266     2 |    0.11    0.01   0.04     0.24   0.01 |     0.41 s
MC footing collapse coarse    616    2602   287 |    0.75    0.77   0.52     3.61   0.08 |     5.73 s
MC footing collapse fine     1244    5178   381 |    2.01    2.46   1.51    10.47   0.14 |    16.59 s
HS footing service coarse     336    1454   584 |    3.11   12.80   0.52     3.35   0.02 |    19.81 s
HS footing service fine       689    2892   563 |    6.27   29.09   1.04     7.28   0.08 |    43.74 s
```

Teşhis (ölçülmüş):

1. **PARDISO her Newton iterasyonunda sıfırdan kuruluyordu** (yerel nesne → sembolik
   analiz + faktörizasyon + yıkım her çağrıda). MC çöküşte sürenin %63'ü, LE'de %59–80'i.
2. **Line-search rezidüel montajı HS'te baskın** (%65–67): her deneme tam bünye
   entegrasyonu (substepping) yapar; HS'in continuum teğeti → lineer global yakınsama →
   adım başına ~15 iterasyon × ~5 line-search denemesi.
3. **COO→CSR sort her iterasyonda** (~%9, MC).
4. Montaj serisel; MKL sequential.

Dürüst not: HS footing GUI yolunda ince mesh'te load_factor=0.844'te kalıyor (çekirdek
test yapısal mesh'te 1.0). Bu, bilinen HS continuum-teğet yakınsama tavanıdır (linchpin);
B4'ün konusu.

## B1 — PARDISO kalıcılık + sembolik faktörizasyon yeniden-kullanımı

`PardisoSolver::factorize` desen-farkındalı yapıldı: aynı çözücü nesnesine sparsity
deseni (row_ptr+col_indices) birebir aynı matris gelirse sembolik analiz atlanır, yalnız
sayısal çarpanlama (phase 22) koşar; desen değişirse (yeni mesh / faz aktivasyonu) tam
analize otomatik düşülür. `build_problem` LinearSolve lambda'ları çözücüyü `shared_ptr`
ile iterasyonlar boyunca yaşatır. SPD/simetrikte sonuç tam yolla birebir; nonsimetrikte
analiz-fazı ölçekleme/eşleme (iparm 11/13) eski değerlerle yeniden kullanılır → round-off
düzeyinde fark olabilir, pivot bozulursa PARDISO hatası tam-analize düşürür (güvenli).

```
LE 40x20 gravity fine       11926   48266     2 |  pardiso 0.24   total  0.40 s   (sembolik 1 kez; tek çözümde kazanç sınırlı)
MC footing collapse coarse    616    2602   287 |  pardiso 1.18   total  3.26 s   (5.73 → 3.26, 1.76×)
MC footing collapse fine     1244    5178   381 |  pardiso 3.55   total  9.53 s   (16.59 → 9.53, 1.74×)
HS footing service fine       689    2892   563 |  pardiso 2.51   total 39.42 s   (pardiso 7.28 → 2.51, 2.9×)
```

load_factor'lar birebir korunur (0.556/0.544/0.976/0.844); HS iterasyon sayısı ±1
(nonsimetrik ölçekleme yeniden-kullanımı round-off'u — yakınsanan sonuçlar değişmedi).

## B2 — CSR desen yeniden-kullanımı (sort'suz değer toplama)

`SparseMatrixBuilder::build_cached(CsrPatternCache&)`: COO girdi imzası (satır,sütun
dizileri) öncekiyle birebir aynıysa sıralama/tekilleştirme atlanır; değerler **build()'in
toplama sırasının birebir aynısıyla** (tarama-sırası haritası) toplanır → sonuç bit-birebir
(floating-point toplama sırası dahi aynı). İmza uyuşmazsa tam build + harita tazeleme.
Newton döngüsünde builder + cache iterasyonlar arası yaşar. Birim test: test_sparse_matrix
`test_build_cached` (bit-birebir + desen-değişimi düşüşü).

```
MC footing collapse fine     1244    5178   381 |  csr 1.47 → 0.22 (6.7×)   total 8.57 s
HS footing service fine       689    2892   563 |  csr 1.04 → 0.15           total 38.79 s
```

İterasyon + line-search sayıları B1 ile BİREBİR aynı (bit-birebir tasarım doğrulandı).

## B3 — Paralel eleman montajı (kalıcı UYUYAN thread havuzu)

İki-fazlı montaj: (1) eleman-başına fe/ke hesabı PARALEL — her eleman yalnız kendi
tamponuna ve kendi Gauss trial durumuna yazar (paylaşılan FP toplama yok); (2) scatter
SIRALI — f_int ve COO düzeni seri yolla birebir aynı ⇒ **sonuç iş parçacığı sayısından
bağımsız bit-deterministik** (ölçüldü: 1/8 thread'de iterasyon sayıları özdeş).
Eleman koordinat+DOF eşlemesi solve başına bir kez toplanır.

**ÖNEMLİ DERS — OpenMP (vcomp) REDDEDİLDİ, ölçümle:** `/openmp` işçileri paralel bölge
sonrası spin-wait yapıyor ve hemen ardından koşan TEK-iş-parçacıklı PARDISO ile çekirdek
yarıştırıyor → PARDISO 3.9→9.25 s (2.4× YAVAŞLAMA), toplam geriledi; OMP_NUM_THREADS=1
ile aynen geri döndü (kanıt). Çözüm: kendi kalıcı havuzumuz (`katai/math/thread_pool.hpp`),
işçiler condition_variable üzerinde UYUR → seri bölümlere sıfır maliyet, harici
runtime/DLL bağımlılığı yok. n<64'te seri (küçük testler birebir eski yol).

```
scenario                    elems    dofs iters | tan-asm res-asm    csr  pardiso  other |    total
LE 40x20 gravity fine       11926   48266     2 |    0.07    0.00   0.05     0.23   0.00 |     0.35 s
MC footing collapse coarse    616    2602   287 |    0.36    0.48   0.06     1.15   0.02 |     2.07 s
MC footing collapse fine     1244    5178   381 |    0.92    1.20   0.19     3.46   0.09 |     5.86 s
HS footing service coarse     336    1454   585 |    1.18    4.33   0.06     1.18   0.01 |     6.76 s
HS footing service fine       689    2892   563 |    2.36    9.74   0.14     2.61   0.02 |    14.87 s
```

## Kümülatif sonuç (B1+B2+B3, i7-7700HQ 4C/8T)

| Senaryo | Taban | B1+B2+B3 | Hızlanma |
|---|---|---|---|
| MC çöküş coarse | 5.73 s | 2.07 s | **2.8×** |
| MC çöküş fine | 16.59 s | 5.86 s | **2.8×** |
| HS footing coarse | 19.81 s | 6.76 s | **2.9×** |
| HS footing fine | 43.74 s | 14.87 s | **2.9×** |

Doğruluk: 85/85 ctest yeşil; tüm load_factor/iterasyon sayıları taban ile birebir.

## B4 — HS sayısal consistent teğet + hibrit politika (LİNCHPIN, 2026-06-13)

Sayısal consistent teğet (Pérez-Foguet & Rodríguez-Ferran & Huerta, CMAME 2000/2001):
tam substep'li güncellemenin Δε'ye göre ileri-fark türevi; pertürbe koşular taban koşunun
alt-adım sayısına SABİTLENİR (nsub_fixed — kritik; ayrıntı hardening-soil-formulation.md §9).
Solver `TangentMode {kNone, kContinuum, kConsistent}` ile ister: line-search rezidüel yolu
kNone (HS'te FD orada hiç koşmaz — "teğet-atlama" küçük kazanımı da bu adımda geldi).

A/B ölçümü (aynı build, geçici KATAI_HS_CONTINUUM anahtarıyla):

| HS GUI senaryosu | continuum | saf consistent | **hibrit (final)** |
|---|---|---|---|
| coarse | 7.5 s, lf=0.976 ✗ | 14.5 s, lf=1.000 | **10.1 s, lf=1.000** |
| fine | 16.2 s, lf=0.844 ✗ | 41.8 s, lf=1.000 | **16.1 s, lf=1.000** |

**HİBRİT POLİTİKA:** artım continuum ile başlar; yakınsamazsa dlam yarılanmadan önce aynı
artım consistent ile yeniden denenir; commit'te continuum'a dönülür. Kalıcı-yükselme
varyantı ölçülüp reddedildi (geçici takılma solve'un kalanını pahalı modda bırakıyor:
structured footing 829 vs 594 iter). Yalnız HS varken devrede → MC/LE iterasyon sayıları
BİREBİR korunur (287/1806, 381/2970 — doğrulandı).

**ANA SONUÇ: GUI-yolu HS yakınsama TAVANI ÇÖZÜLDÜ** (load_factor 0.976/0.844 → 1.000/1.000,
neredeyse continuum maliyetine). Dürüst not: consistent teğet düz rejimde kuadratik hız
BEKLENTİSİNİ vermedi (tol=1e-2'ye continuum da 3-4 iterde ulaşıyor; aktif-set/spektral
kıvrımlar FD kolonlarını geçişlerde harmanlıyor) — değeri zor-rejim robustluğunda.

Kümülatif oturum kazanımı (taban → B4): MC çöküş fine 16.6→6.1 s (2.7×, bit-birebir);
HS footing fine 43.7 s @ lf=0.844 → 16.1 s @ **lf=1.000** (2.7× + tam yakınsama).

## Sıradaki adımlar

- Undrained yüksek-yük tavanını (test_hs_undrained_bvp kısıdı) hibritle yeniden ölç —
  aynı kökten çözülmüş olmaya aday.
- MKL threaded (intel_thread + libiomp5) yalnız faktörizasyon yeniden baskınlaşırsa
  (büyük mesh'lerde). Stress recovery / K0 montajı gibi tek-seferlik yollar ikincil.
- GUI'de hesap-süresi/iterasyon raporu (SolveResult.timings hazır).


## Nonlineer dinamik maliyet çarpanı (2026-07-20, ölçülü)

Aynı model (SH kolonu, MC, 600 Newmark adımı, test_dynamic_gui): **lineer 0.09 s, nonlineer 6.07 s
→ 70×.** Sebep yapısal: lineer yol K_eff'i BİR kez faktörler; nonlineer yol her Newton
iterasyonunda K_T'yi yeniden kurar+faktörler (nonsimetrik mtype=11). Bu çarpan ürünün kendi
yardımında beyan edilir ("much slower, opt-in"); modifiye-Newton (teğeti adımda dondurma)
seçeneği gelecekte bu çarpanı düşürmek için kayıtlı adaydır (ROADMAP Track 1 kuyruğu).

---

# Build-time baseline (roadmap section 5.2)

Distinct from the runtime campaign above: these rows measure COMPILE cost, because a
one-hour clean build is a correctness risk (it discourages the full-suite runs that
catch silent-wrong). Method: per-translation-unit wall times from Ninja's own build
log (`build/msvc-rwdi/.ninja_log`, latest entry per output) — every build measures
itself for free. Machine: 8 hardware threads, MSVC RelWithDebInfo, shared precompiled
header, dynamic MKL, ccache; effective parallelism ~6–8 under the JOB_POOLS caps.

## Batch 1 — the driver bodies leave the headers (2026-08-01)

Measured before: the driver family (`katai/jobs`) carried its function bodies
header-inline — `solve_gravity_le` alone ~1,200 lines — and every one of the ~26
driver-including test TUs re-generated code for those bodies and the engine machinery
they pull in, at ~205–260 s per TU. The bodies (1,636 lines) moved to
`kernel/jobs/src/{driver,mesh_builder,flow_driver}.cpp`, `katai_jobs` became a
compiled static library, and the headers keep the types, small helpers and
declarations.

| Quantity | Before | After | Change |
|---|---|---|---|
| Sum of TU compile times, whole tree | 9,741 s | 5,547 s | **−43%** |
| Sum over `tests/` (122 TUs) | 8,689 s | 4,367 s | −50% |
| Worst driver-test TU (`test_dynamic_gui`) | 261.6 s | 69.7 s | −73% |
| Typical driver-test TU | ~205–235 s | ~55–70 s | ~3–4× |
| The driver body itself (`src/driver.cpp`, paid once) | — | 113.2 s | ×1, no longer ×26 |
| Estimated clean-build wall clock at ~7× parallelism | ~23 min | ~13 min | −10 min |

Correctness: full rebuild against the split headers, suite 126/126; the architecture
gate stays green (the `src/` bodies are module-private).

**Next targets, measured (top offenders after batch 1):** `nonlinear_solver.cpp`
~120 s, `dynamics_nonlinear.cpp` ~119 s, `test_ssi_dynamics` 110 s,
`test_consolidation_plastic` 102 s, `test_coupled_flow_plastic` 97 s — these point at
`katai/analysis`'s template-heavy headers (nonlinear solver, dynamics), the §5.2
item-2 batch over the engine modules. The <10-minute clean-build target needs that
batch; this one removed the suite's widest multiplier.

## Batch 2 — the driver headers carry the declaration closure only (2026-08-02)

The ~20 engine headers the solve bodies consume (phase strategies, dynamics,
consolidation, interfaces, assembly, solvers) moved from `driver.hpp` /
`mesh_builder.hpp` / `flow_driver.hpp` into the `src/*.cpp` files, together with every
body-only helper (measured first: zero external consumers for each). The headers now
hold types, enums, four externally-used helpers and declarations. Fallout was three
files: two tests gained the direct engine include their direct engine use always
deserved, and the facade gained `backend_name()` (result provenance — the About box
already showed it, and `katai info` will want it).

**Measured honestly, two different outcomes:**

| Quantity | Result |
|---|---|
| Clean-build TU-time sum | ~neutral (5.5–5.7 ks; the shared PCH is Eigen+std only, and raw header PARSE was never the fat cost — batch 1's codegen was) |
| **Incremental cone: touch `phase_solver/dynamic.hpp`** | **rebuilds exactly 1 TU (`src/driver.cpp`, ~2 min) — before the diet it rebuilt every driver consumer, ~35 TUs ≈ half an hour** (measured with a touch + `ninja -n` dry run) |

The incremental cone is the §5.2 target that matters day to day ("incremental build
after one `.cpp` edit under 30 seconds" — an engine-header edit is now one library TU
away from it), and it is what the next batch inherits: engine-header work stops
costing a suite-wide rebuild. Suite 127/127; architecture gates green.