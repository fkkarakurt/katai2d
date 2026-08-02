# KATAI 2D — Doğrulama Özeti (PLAXIS 2D + analitik kıyaslar)

Bu tablo KATAI'nin her özelliğini, sektör programlarının da doğrulandığı **analitik kapalı-form**
çözümlerle ve/veya **PLAXIS 2D** sonuçlarıyla karşılaştırır. **Mantık:** analitik referansı olan bir
problemde PLAXIS de o analitiği tutturur; dolayısıyla KATAI'nin analitiğe <%1-2 yakınlığı, PLAXIS ile
**birebir/çok-yakın** olduğu anlamına gelir. Her satır otomatik bir regresyon testidir (`ctest`; current
baseline **120/120 green** — see the §0 baseline stamp).

Kaynaklar: PLAXIS 2D 2025.1 Scientific/Reference/Material Models Manual; klasik geoteknik analitikleri
(Prandtl, Reissner, Lamé, Terzaghi, Flamant, Boussinesq, Rankine, Thiem, Harr). Detaylar: aşağıdaki testler
ve `docs/references/*-formulation.md`.

## 0. Doğrulama & Onaylama (V&V) çerçevesi

Bu belge KATAI 2D'nin **resmî V&V kaydıdır** ve NAFEMS/AGS *"Validation and Use of Geotechnical
Software"* ilkeleriyle uyumlu yapılandırılmıştır. İnsan hayatını etkileyen tasarım üreten bir araç
olduğundan doğrulama kapsamı **sürüm kapısıdır** (bkz. `docs/internal/roadmap-v0.3-mvp.md` §7).

**Doğrulama yöntemi (üç katman):**
1. **Analitik onaylama** — kapalı-form çözüme sahip problemlerde mutlak hata bandı (< %1–2, çoğu round-off).
2. **Benchmark kıyası** — analitiği olmayan problemlerde PLAXIS 2D / literatür FE sonuçlarıyla karşılaştırma.
3. **Regresyon** — her satır bir `ctest` testidir; her build'de otomatik yeniden koşar (sessiz sapma yakalanır).

**Risk sınıflandırması (NAFEMS QSS — Vital / Important / Advisory):** özelliğin çıktısının bir tasarım
kararında taşıdığı risk düzeyi. *Vital* = doğrudan can güvenliği/taşıma gücü/stabilite; *Important* =
tasarımı etkileyen ama dolaylı; *Advisory* = yardımcı/gösterim.

| Bölüm | Kapsam | Risk sınıfı |
|---|---|---|
| §1 Plastisite limit yükleri (taşıma gücü, itki, şev FoS) | Çökme/stabilite | **Vital** |
| §2 Lineer elastik analitik (oturma, gerilme dağılımı) | Servis-yük tepkisi | Important |
| §3 Su — sızma + konsolidasyon (uplift, heave, oturma-zaman) | Kaldırma/heyelan/konsolidasyon | **Vital** |
| §4 Yapısal elemanlar + kuplaj (perde, ankraj, kazık) | Yapısal güvenlik | **Vital** |
| §5 PLAXIS 2D sonuç-kıyasları (FE-vs-FE) | Doğrulama çapraz-kontrol | Important |
| §6 GUI compute-path (kullanıcının gerçek "Calculate" yolu) | Uçtan-uca doğruluk | **Vital** |

**Current baseline stamp (2026-07-30):** `120/120 tests passed, 0 failed` on the default (PARDISO)
build (160 s) · `118/118` on the `portable` Eigen build (242 s; the two PARDISO adapter unit tests
exist only where the adapter does) · `90/90` on the `engine` composition (no project schema, no GUI;
the 30 schema-dependent tests are skipped *by name* at configure time) · exit code 0 everywhere ·
`scripts/check_composition.ps1` walks all three rows.
Environment: Win11, MSVC RelWithDebInfo (Ninja), Intel oneMKL 2026.0 sequential dynamic (PARDISO), ccache.

**Kapsam & henüz-kapsanmayan:** aşağıdaki §1–6 mevcut çekirdeği kapsar. v0.3'te eklenecek özellikler
(EC7/TBDY 2018 tasarım katsayıları, dinamik/sismik, Soft Soil/Cam-Clay/Hoek-Brown bünye) **kendi V&V
satırları eklenmeden "biter" sayılmaz** — DoD sözleşmesi: `docs/internal/roadmap-v0.3-mvp.md` §6.

## 1. Plastisite limit yükleri (en keskin test)
| Problem | KATAI | Referans (analitik) | Hata | Test |
|---|---|---|---|---|
| Prandtl şerit temel Nc (φ=0) | 5.20 (tri15) | 2+π = 5.14 | %1.1 | test_tri15_prandtl |
| Reissner Nc (φ=20°) | 15.13 | 14.84 | %2.0 | test_bearing_phi |
| Silindir çökme (axisym) | — | 2c·ln(b/a) | %0.3 | test_axisym_collapse |
| Rankine aktif itki | 300.0 | ½Ka γH² = 300 | round-off | test_earth_pressure |
| Rankine pasif itki | 2700 | ½Kp γH² = 2700 | %~0 | test_earth_pressure |
| Rankine at-rest (K0) | 450.0 | ½K0 γH² = 450 | round-off | test_earth_pressure |
| Şev FoS (φ-c reduction/SRM) | 1.01 | ~0.99 (Bishop) | %2.0 | test_slope |
| **MC çekme kesmesi (Rankine kapağı)**: 6 bölge tersinmesi + kapalı formlar (tek-eksenli plato, T-apeks, TT-kenarı) + 8000'lik uygunluk/idempotens taraması + 4×4 J vs merkezî-FD | round-off / <1e-4 (J) | MMM Denk. 3-11 + formülasyon §7 kapalı formları | kesin | test_mohr_coulomb |

## 2. Lineer elastik analitik (kapalı-form)
| Problem | KATAI hata | Referans | Test |
|---|---|---|---|
| Lamé kalın silindir (axisym σr,σθ) | max 1e-5 | Lamé | test_axisym_cylinder |
| Boussinesq şerit yük σz | <%5 | Boussinesq | test_boussinesq |
| Flamant nokta/çizgi yük σz=2P/πz | <%6 (kenar) | Flamant | test_flamant |
| Konsol kiriş (plate) PL³/3EI+PL/kGA' | round-off (tek tri15-eleman) | Timoshenko | test_plate, test_plate5 |

## 3. Su — sızma + konsolidasyon (PLAXIS Sci.Man §3,§4)
| Problem | KATAI hata | Referans | Test |
|---|---|---|---|
| 1D Darcy kolonu (head, debi) | round-off | analitik | test_seepage |
| Radyal akış (Thiem) head/debi | order O(h³) / 1e-6 | Thiem | test_seepage |
| Baraj-altı uplift h(x)=(H/π)arccos | %0.49 | Harr arccos | test_seepage |
| Dupuit/Charny iki-rezervuar debi | %0.08 | Charny | test_seepage |
| **Terzaghi 1D konsolidasyon U(Tv)** | **<%2** | Terzaghi serisi | test_consolidation |
| **Transient akış izokron u(Z,Tv)** (W1) | **~%1** | Terzaghi Fourier serisi | test_transient_flow |
| **Transient ani-head adımı h(d,t)** (W1) | **<%0.4** | Carslaw-Jaeger erfc | test_transient_flow |
| **Transient → kararlı-hal (steady seepage)** (W1) | **9e-14** | lineer head profili | test_transient_flow |
| **van Genuchten/Mualem analitik kapasite=FD** (W2a) | **~1e-11** | van Genuchten (1980)/Mualem (1976) | test_water_retention |
| **Richards infiltrasyon kütle-korunumu** (W2b) | **1.000000** | Celia (1990) kütle dengesi | test_unsaturated_flow |
| **Richards Philip √t sorptivite** (W2b) | I(8)/I(2)=2.20 | Philip √t (lineer=4.0) | test_unsaturated_flow |
| **Fully-coupled (W3) doygun = konsolidasyon** | **max\|Δ\|=0 (kesin)** | solve_consolidation (Biot) | test_coupled_flow |
| **Fully-coupled (W3) doymamış Bishop aktivasyon** | S_eff=0.8575 birebir | van Genuchten + resaturasyon | test_coupled_flow |
| **GUI FullyCoupled fazı = Terzaghi** (W4) | <%5 | Terzaghi U(Tv) | test_water_gui |
| **GUI TransientFlow fazı** (W4) | akış+resaturasyon | basınçlandırma→pore | test_water_gui |
| **HSsmall** degradasyon (G_s,G_t,cut-off) | **round-off** | Hardin-Drnevich (MMM §7) | test_hssmall |
| HSsmall FE (E0 küçük-şekil. → HS degrade) | exact gate | MMM §7 | test_hssmall |
| MMS yakınsama (V&V) | tri6 O(h³), tri15 O(h⁵) | Roache | test_seepage |

## 4. Yapısal elemanlar + kuplaj (PLAXIS Sci.Man §7, Ref §5.6)
| Problem | KATAI hata | Referans | Test |
|---|---|---|---|
| Ankraj (fixed-end + node-node) elastoplastik | round-off | kapalı-form | test_anchor, test_anchor_plastic |
| Geogrid tension-only + N_p | round-off | kapalı-form | test_geogrid |
| Interface Coulomb çökme λ* | %0.25–1 | limit-analiz | test_interface |
| Interface K0 install (wished-in-place) | u=0 (4e-12) | öz-denge | test_wall_k0_excavation |
| Ters izoparametrik harita (eğri eleman) | 1.4e-15 | round-trip | test_point_location |
| **Embedded beam eksenel (sabit zemin)** | **round-off** | (P/EAλ)coth(λL) | test_embedded_beam |
| **Embedded beam yük transferi (N_s)** | **round-off** | Σ skin = P | test_embedded_beam |
| **Embedded beam nihai kapasite** | **−%0.1** | Q_skin+Q_base | test_embedded_beam |
| **Plate M-N mafsalı dönüş haritası (MMM §18.3 elması)** | **kesin** (f=0.00, CPP oracle birebir) | brute-force enerji-normu izdüşümü | test_plate_plastic |
| **Plate mafsalı BVP: onset 2·Mp/s_g + tam-Mp doyması** | onset braketi ±%10; doyma **tam** (100.0000000000) | kapalı-form (Gauss-kollu) | test_plate_plastic |
| **Plate M-N elmas yüzeyi (eksenel+eğilme BVP)** | f = 0.000e+00; plato ~Mp(1−N/Np) | elmas kesiti + Np-sınırsız kontrol | test_plate_plastic |

## 5. PLAXIS 2D sonuç-kıyasları (FE-vs-FE)
| Problem | KATAI vs PLAXIS | Not | Doküman |
|---|---|---|---|
| Kantilever palplanş kilde, iyi-gömülü | −%7 … +16% | tri15, undrained-B | wall-benchmark-plaxis.md |
| (derin-gömülü) | domain yakınsamasıyla PLAXIS'e | sınır≥kırılma kaması | wall-benchmark-plaxis.md |
| Embedded beam kapasite | by-construction tutarlı | T_max/F_max = PLAXIS girdisi | — |

## 6. GUI compute-path (entegre `build_problem`) doğrulaması — `gui-pipeline-validation.md`
Yukarıdaki 1–5 çekirdek kernel'leri *doğrudan* doğrular. Bu satırlar kullanıcının GUI'de bastığı
**"Calculate"in tam yolunu** (mesh_from_project → solve_gravity_le → nodal recovery) niceliksel sınar.
| Problem | KATAI (GUI) | Referans | Hata | Test |
|---|---|---|---|---|
| 1B ödometre (öz-ağırlık + surcharge) | FE | `−γH²/2E_oed` / `−qH/E_oed` | <%2 | test_solve |
| K0 procedure (max\|u\|, efektif σ') | 0.0 / FE | jeostatik denge | 1e-13 | test_gui_solve, diag_k0 |
| MC Prandtl `N_c` (φ=0) GUI yolundan | 4.84 | 5.14 | −%5.9 (konservatif) | study_gui_validation |
| MC Reissner `N_c` (φ=20°) GUI yolundan | 13.94 | 14.84 | −%6.1 | study_gui_validation |
| MC = LE elastik rejimde | birebir | LE ödometre | <0.1% | study_gui_validation |
| **Axisym K0 silindir** (max\|u\|, σ_zz, σ_rr) | 0.0 / FE | jeostatik | 1e-13 | test_axisym_gui |
| **Axisym gravity ödometre** `−γH²/2E_oed` | −0.06686 | −0.06686 | <0.01% | test_axisym_gui |
| **Axisym radyal genişleme** `u_r(R)=νqR/E` | birebir | hoop ε_θ=u_r/r | round-off | test_axisym_gui |
| **MC çekme kesmesi tam yolda** (5 m düşey kazı, c=25 kil; MMM §3.3.10 hendek senaryosu): kret arkası max σ₁ | AÇIK 0.40 kPa / KAPALI 7.71 kPa | σ_t=0 kapağı (recovery payı) | ikili tanık | test_input_audit |
| **Yapısal öz-ağırlık statikte** (plate w, kazık γA/Ls): tutarlı düğüm kuvvetleri el-integralleri (Simpson L/6-L/6-2L/3; Boole 7-32-12-32-7/90; eğik toplam −wL; trans_dof; kazık kendi DOF'ları) + K0'da w=0 kimliği / w=8.3 oturma tanığı / ağırlıklı nil-faz no-op | kesin (1e-11) / 0.000e+00 / 1.80 mm / 0.000e+00 | el integralleri + f/B tek-sayım sözleşmesi | kesin | test_input_audit |
| **Soft Soil Creep GUI yolu** (Aşama 3): K0 tohumu 0.000e+00; NC kolonda zaman sınıfları (T=0 κ*-sınıfı; T=τ ~λ* izokronu %25 sınıf bandı); **u(100g)−u(1g) = μ*·H·ln100 %+0.1** (SumMstage zaman-paylaştırma "gençlik açığı" farkta düşer); SSC+Undrained açık red; IO birebir | +%0.1 (creep kuyruğu) | MMM §11 izokron kapalı formları | kesin sınıf | test_soft_soil_gui |
| **NonPorous drenajı** (tam batık kolon, su tablası yüzeyde): σ_yy(orta) = −γ_unsat·z **birebir −96.000** (kaldırmasız, pore'suz TOPLAM) + K0 kimliği **0.00e+00** (tohum+gravite+pore-dışlama tutarlı); Drenajlı ikiz −60.760 = −γ'·z tanığı; konsolidasyonda açık red | kesin | PLAXIS gözeneksiz kuralı (γ_unsat, pore yok, geçirimsiz) | kesin | test_input_audit |
| **SSC tutma deseni** (zamanlı yük + zamanlı NİL zinciri): yaşlanmış zincir ideal 100-gün izokronu λ*I+μ*Hln100 üzerine **+%0.7**; rampalı koşu daha genç → daha az oturur (sıralama pinli). **SSC×Konsolidasyon DOĞRULANAMADI → AÇIK RED** (oturma zaman-adımına bağımlı: 50/100/200 adım −0.519/−0.506/ıraksama; kuyruk −%37) — Biot×creep araştırması açık iş kalemi | +%0.7 / red pinli | izokron kapalı formu + adım-bağımsızlık ölçümü | kesin | test_soft_soil_gui |
| **Sistematik matris** K0 admissibility 20/20 | max\|u\|=0 | LE/MC×tri6/tri15×5 yapı | round-off | test_gui_matrix |
| **Sistematik** her yapı hesaba katılıyor + tri6≈tri15 + su buoyancy | FE | değişmezler | — | test_gui_matrix |
| **Safety (φ-c reduction) GUI** şev FoS | 1.009 | Bishop ~0.99 | %1.9 | test_safety_gui |
| Mesh domain bütünlüğü (yük/yapı taşma yok) | 0 taşma | CDT | exact | test_mesh_domain |

**Bulgular:** (a) K0 procedure DOĞRU hesaplıyor; "GUI'de çalışmıyor" yalnız gösterim hatasıydı (varsayılan
alan + deplasman ölçeği), düzeltildi (commit 74057ec). (b) GUI çöküş yükünü ~%6 alttan (konservatif) tahmin
eder — artımsal limit-analizin yük-adımı granülaritesi (`min_dlam`); elastik tepki TAM. Detay + kalan denetim
listesi: `gui-pipeline-validation.md`.

## 7. Tasarım kodu uyumu — Eurocode 7 (EN 1997-1) [Risk sınıfı: **Vital**]

EC7 **malzeme-katsayılı** yaklaşımlar (DA1-C2, DA3) mevcut φ-c reduction makinesiyle çözülür ve
**PLAXIS "Design Approaches" ile ALGORİTMA-ÖZDEŞ** (Bentley/PLAXIS: cohesion/friction/dilatancy
partial-factor ile azaltılır, yük partial-factor ile artırılır). Değerler `docs/references/design-codes-ec7-tbdy.md`'de kilitli.

| Problem | KATAI | Referans (analitik/resmî) | Hata | Test/Study |
|---|---|---|---|---|

## 8. Kusur kaydı

**✅ RESOLVED — interface 2026-07-28, and every direct caller moved onto it the same day (risk
class: Vital): a direct solver can answer a singular system instead of refusing it.** A rank-deficient
nonsymmetric matrix does not make PARDISO fail. It perturbs the tiny pivots, completes the
factorization with `error == 0`, and the following solve returns a finite solution that satisfies
nothing. Since a singular global stiffness is what an insufficiently restrained model produces, the
failure mode is a plausible displacement field where an error is required — and it is invisible,
because the driver's "insufficiently restrained" message only appears when the solver throws.
**Measurement:** a deliberately rank-deficient 4×4 system (two identical rows, dense rank 3 < 4)
returned a finite answer with no error code, while Eigen's sparse LU refused the same matrix outright.

**The first attempted fix was wrong, and measuring it is what showed why.** Refusing whenever
`iparm(14)` reports a perturbed pivot broke three previously passing tests — `test_embedded_beam`,
`test_axisym_collapse`, `test_earth_pressure` — all of which produce correct, validated results while
PARDISO perturbs pivots as routine numerical robustness. **Perturbation count is therefore not a
singularity signal**, and it is now exposed as a diagnostic only (`PardisoSolver::perturbed_pivots`).

**Fix:** acceptability is judged from the solution rather than from the factorization. Every solve
through `katai::linsolve::DirectSolver` is verified against the system it claims to solve: relative
infinity-norm residual ‖Ax−b‖/max(‖b‖,1), refused above 1e-6. The check lives in the non-virtual part
of the interface, so no backend can be written that skips it, and it is backend-independent by
construction. Measured: the singular 4×4 yields residual 0.20 and is refused; healthy systems report
below 1e-10. Cost is one sparse matrix-vector product per solve against a factorization orders of
magnitude more expensive.

**Closed for the product path (2026-07-28):** every direct call site — the staged-analysis driver
`build_problem.hpp`, the flow driver `build_flow.hpp`, the GUI's strip-load demo, and every test that
solved through `math::PardisoSolver` — now solves through `katai::linsolve::DirectSolver`, so the
residual check guards the product's real solve path, not only the interface's own unit test. No C++
source in the tree names a backend or conditions on `KATAI_WITH_MKL`. Alongside, a refused solve
became a vocabulary: `SingularSystem` (recoverable — the tangent is rank-deficient along a collapse
mechanism) is distinguished from `SolveError` (a fault), `solve_nonlinear` catches only the former,
counts it in `NewtonResult::refused_solves`, and abandons the increment so the outer loop cuts back —
the collapse-load path now works by contract where it previously worked by accident (the meaningless
finite vector happened to fail the line search). Measured: the full suite passes on both backends —
119/119 on the PARDISO build, 117/117 on the `portable` Eigen build (the two PARDISO adapter unit
tests exist only where the adapter does). The architectural remainder — the driver still *living* in
the application layer — is Stage B work and is tracked in the architecture exceptions, not here.
**Permanent tripwire:** `test_linsolve` asserts that *both* backends refuse the singular system, and
the fixture verifies its own rank deficiency with a dense factorization, so it cannot decay into a
vacuous test.

**✅ RESOLVED (measured and fixed 2026-07-28): the documented aliasing contract of
`PardisoSolver::solve` was false (risk class: Important).** The header states that `rhs` and `solution`
may point at the same array. PARDISO with `iparm[5] = 0` writes the solution while still reading the
right-hand side, so an aliased call silently corrupts the result. **Measurement:** on a 20-DOF SPD
system an in-place solve deviated from the out-of-place solve by max|Δx| ≈ 62 — a wrong answer with no
diagnostic. **Fix:** the wrapper takes a scratch copy of the right-hand side, and only when the two
pointers actually overlap, so the common non-aliased path is unchanged. **Permanent tripwire:**
`test_linsolve` performs an in-place solve on every backend and compares it with the out-of-place
result. Both defects were found by cross-checking two independent solver backends against a dense
factorization — neither is visible from a single backend's own output.

**Recorded 2026-07-30 — the Stage A2 split batch changed no physics and moved no number.** `core`
became `katai/materials` + `katai/fem` + `katai/analysis`, and the driver's inline SoilModel switch
became the constitutive registry (`katai/materials/registry.hpp`): same construction arithmetic,
reached through one named seam. The claim that this is behaviour-preserving is *measured*, not
asserted: every validated pin in the suite passed unchanged on both backends after the restructuring
(120/120 PARDISO, 118/118 Eigen; the +1 over the previous stamp is `test_material_registry` itself,
which pins build == the former inline construction field-exact and the refusal messages verbatim).
The `engine` composition (90/90 with the schema and GUI absent) additionally proves the physics
tests do not depend on anything above the engine.

**✅ ÇÖZÜLDÜ (2026-07-19 ölçüldü → aynı gün düzeltildi): statik faz zinciri yapısal durumu taşımıyordu
(Risk sınıfı: Vital).** `solve_phases` fazlar arası yalnız committed Gauss gerilmesini taşıyordu;
yapısal elemanlar her zincirli fazda u=0 + sıfır plastik hafızayla başlıyor, SumMstage dengesizliği
ebeveynin yapısal traksiyonlarını YENİDEN ramp'lıyordu. **Ölçüm:** değişikliksiz nil fazda duvar
momenti %32 kayıyordu (max|ΔM|=1.95 / |M|=6.03 kNm/m) ve K0→Plastic el-değişimi K0'ın kendi yakınsamış
duvar durumundan ~5× sapıyordu. **Düzeltme:** Track 1a deseni statiğe genellendi — `solve_nonlinear`
opsiyonel `StructuralInit` (ebeveyn deplasman datumu + committed plastik durum) alır; baseline B tam
yapısal iç kuvveti (f_s0, paylaşılan montajcıyla) içerir → residual(0)=0 ebeveynin KENDİ dengesiyle
sağlanır, ramp yalnız gerçek konfigürasyon değişikliğini taşır. **Kimlik (test_staged_struct_carry,
artık kalıcı tripwire, 1e-6 bağıl tolerans):** nil faz |u|=1.3e-11 m, duvar ΔM=2.2e-8 kNm/m, arayüz
Δτ=6.8e-7 kPa — gerçek no-op. Sıfır-yapılı zincirler bit-birebir korunur. Uçtan uca kimlik
(K0 → statik taşıma → sürşarj → dinamik taşıma → sürücüsüz sarsıntı) istasyon-istasyon 0.000e+00
(test_dynamic_gui). Plastik-olay dedektörü pikometre eşikli (1e-12 m): tam akma sınırında
yeniden-değerlendirmenin ulp-gürültüsü [SLIPPING] rozetini tetikleyemez.
| DA set tabloları (A/M/R × DA1-C1/C2/DA2/DA3) | birebir | EN 1997-1 Annex A önerilen | exact | test_design_code |
| Malzeme-katsayısı (c/γ_c', tanφ/γ_φ', Tresca c_u/γ_cu) | birebir | EN 1997-1 §2.4.7.3 | round-off | test_design_code |
| **DA3 şev ODF = FoS_karakteristik / γ_M** | 0.8074 | 1.0091/1.25 = 0.8073 (kesin özdeşlik, mesh-bağımsız) | **%0.017** | study_design_ec7 |
| DA3 şev FoS_karakteristik | 1.009 | Griffiths & Lane / Bishop ~0.99 | %2 | study_design_ec7 |
| EC7 DA3 tasarım hükmü (ODF<1 ⇒ güvensiz) | doğru | ODF≥1.0 kriteri (EN 1997-1) | — | study_design_ec7 |

**Not (TBDY 2018):** direnç-katsayılı `E_d ≤ R_d = R_k/γ_Rv` yaklaşımı rapor/kontrol katmanında;
γ_Rv değerleri birincil standarttan (TBDY 2018 Tablo 16.1) pinlenince V&V satırı eklenecek.

## Değerlendirme
- **Analitik-destekli özellikler (1–4): KATAI <%1-2 (çoğu round-off).** Bu, PLAXIS'in de tutturduğu altın-
  standart analitiklerdir → **PLAXIS ile birebir/çok-yakın.**
- **Doğrudan PLAXIS sayı-kıyası (5):** iyi-koşullu vakalarda mühendislik bandında; sapmalar kurulum
  (domain/mesh) kaynaklı, fizik değil (domain büyüdükçe PLAXIS'e yakınsar — kanıtlandı).
- **Algoritma uyumu:** stress-point return-mapping + tutarlı teğet = PLAXIS/FLAC ile aynı mimari
  (Sci.Man §2 deformasyon teorisi birebir). Eleman ailesi tri6+tri15 (PLAXIS varsayılanı), yapısal 5-düğüm.

**Açık kalan PLAXIS-parite işleri** (`plaxis-gap-analysis.md`): HSsmall, embedded-beam ISF (load-settlement
rijitliği), graded-mesh wall (derin vakalar verimli <%5), convergence-criteria hizalama. **Su hikayesi:**
W1 transient doygun akış ✅ + W2 doymamış van Genuchten/Mualem + Richards ✅ + W3 fully-coupled flow-deformation
(LE iskelet, Bishop χ=S_eff; doygun-limit bit-birebir) ✅ + W4 GUI entegrasyonu (TransientFlow + FullyCoupled faz
tipleri) ✅; kalan W3-takip (elastoplastik iskelet, Liakopoulos benchmark) — bkz `transient-unsaturated-flow-formulation.md`.
