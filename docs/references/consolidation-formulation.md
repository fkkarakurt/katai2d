# Biot Konsolidasyon — Formülasyon (Faz A.2)

Zaman-bağımlı kuplajlı deformasyon–su akışı (Biot 1941). Undrained (anlık) ve drained (nihai) ARASI:
boşluk suyu basıncı zamanla sönümlenir, efektif gerilme artar, oturma gelişir. PLAXIS'in en sık
kontrol edilen analizlerinden; klasik **Terzaghi 1D** kapalı-form doğrulaması (altın standart).

**Kaynaklar (kilitli):** PLAXIS 2D 2025.1 Scientific Manual §4 (Eq 4-1…4-32); Biot (1956); Verruijt;
Terzaghi (1943) 1D konsolidasyon analitiği. Konvansiyon: tension-pozitif (basınç<0), tüm kod tabanıyla aynı.

## 1. Yöneten denklemler (Sci.Man §4.1)
Terzaghi prensibi: σ = σ' + m(p_steady + p_excess), m=(1,1,1,0,0,0)ᵀ. Bünye σ̇'=M ε̇.
Süreklilik (Eq 4-14, p=excess, steady çıkarılmış): ∇ᵀ(k∇p/γw) + mᵀ ∂ε/∂t − (n/Kw) ∂p/∂t = 0.

## 2. FE blok sistemi (Sci.Man §4.2, Galerkin, aynı N hem u hem p)
u=Nv, p=N pₙ, ε=Bv. Eleman matrisleri:
```
K = ∫ Bᵀ M B dV         (rijitlik, 2n×2n)
L = ∫ Bᵀ m N dV         (kuplaj, 2n×n)
H = ∫ (∇N)ᵀ (k/γw) (∇N) dV   (akış/permeabilite, n×n)   [seepage He'nin γw-ölçekli hali]
S = ∫ (n/Kw) Nᵀ N dV    (depolama/sıkışabilirlik, n×n)
```
PLAXIS Kw/n = 3(νu−ν)/((1−2νu)(1+ν)) · K_skeleton, νu default 0.495 (Eq 4-17; undrained wrapper'la AYNI).
∇N satırları strain_displacement B'sinden (seepage'deki gibi: G(0,i)=B(0,2i), G(1,i)=B(1,2i+1)).

## 3. Zaman integrasyonu (Sci.Man §4.2, Eq 4-18…4-20; α=1 tam-implicit)
```
[K    L ] [Δv ]   [0    0  ] [v0 ]   [Δf_n   ]
[Lᵀ  −S*] [Δpₙ] = [0  ΔtH ] [pₙ0] + [Δt qₙ* ]      S* = αΔtH + S = ΔtH + S
```
Satır-1: K Δv + L Δpₙ = Δf_n   (denge; konsolidasyon fazında Δf=0).
Satır-2: Lᵀ Δv − S* Δpₙ = ΔtH·pₙ0 + Δt qₙ   (süreklilik; kapalı sınır qₙ=0 → ΔtH·pₙ0).
Her adımda 2×2 blok çözülür → v+=Δv, p+=Δp. t→∞'da Δp→0, p→0 (tam drenaj).

## 4. Kritik zaman adımı (Sci.Man §4.4, Vermeer&Verruijt 1981)
Δt < eşik altında doğruluk düşer + eşit-mertebeli u–p tri6'da erken-zaman dama-tahtası boşluk basıncı
belirir (bkz `docs/validation/lbb-undrained-checkerboard.md`). Δt_crit = h²γw/(η kᵧ)(1/Eoed + n/Kw),
η=80 (tri15) / 40 (tri6), h = eleman boyu. Eoed = E(1−ν)/((1+ν)(1−2ν)). **Çekirdek:**
`consolidation_critical_dt` (analysis/consolidation.hpp) + `oedometer_modulus`. Çözücü eşit adım
Δt=süre/adım kullanır (Δt'yi ZORLAMAZ); GUI faz-editörü Δt<Δt_crit ise UYARIR
(`consolidation_step_warning`, build_problem.hpp) — oturmalar her koşulda güvenilir, yalnız erken-zaman
pore alanı etkilenir. Kalibrasyon: `test_consolidation` `test_critical_dt` (LBB study'siyle birebir).

## 5. Terzaghi 1D doğrulama (analitik altın standart)
Yanal-konfine kolon (H_dr drenaj yolu), üstte drenajlı (p=0), altta geçirimsiz, ani yük q.
t=0⁺ undrained: p=q üniform, σ'=0, oturma=0. t→∞ drained: p=0, σ'=q, oturma s∞=qH/Eoed.
Konsolidasyon derecesi U(Tv), Tv = cᵥt/H_dr², cᵥ = k Eoed/γw (Kw→∞):
```
U = 1 − Σ_{m=0}^∞ (2/M²) exp(−M² Tv),  M=(2m+1)π/2
pratik: Tv=(π/4)U² (U≤%60); Tv=−0.933 log(1−U)−0.085 (U>%60)
```
FE U(t) = oturma(t)/s∞ VEYA 1−ortalama_pore/q → Terzaghi eğrisine oturmalı (<%2-5).

## 6. KATAI çekirdek uygulaması (TAM — test_consolidation)
- `analysis/consolidation.hpp`: assemble K/L/H/S (eleman-generic tri6/tri15, mevcut B + seepage G yeniden);
  birleşik DOF (öteleme 2/düğüm DofMap + ek pore DOF 1/düğüm); drenaj sınırı = pore DOF sabit (p=0).
- İlk faz LinearElastic (§4.1-4.2). Elastoplastik konsolidasyon (§4.3, residüel iterasyon) sonra.
- Doğrulama `test_consolidation`: Terzaghi 1D U-Tv birkaç Tv'de + nihai oturma qH/Eoed (dense LU, altın standart).

## 7. GUI faz entegrasyonu (TAM — `PhaseType::Consolidation`, test_consolidation_gui)
PLAXIS "Consolidation" hesap fazı: `build_problem.hpp` içinde `InitialPhase::Consolidation` branch'i
(`solve_gravity_le` içinde, Safety'ye paralel). Çalışma mantığı (v1):
- **Yük artımı t=0+'da uygulanır:** fazın konfigürasyon dengesizliği **dF = f − B** (f = aktif gravity+yükler,
  B = committed iç kuvvet = SumMstage; gravity zaten B'de dengeli → dF = YENİ sürşarj/dolgu) çekirdeğin yeni
  opsiyonel `load_increment` parametresiyle İLK zaman adımında uygulanır. S≈0 (sıkışmaz su) → t=0+ tepkisi
  UNDRAINED'dir: yük fazlalık boşluk basıncı üretir (Skempton B≈1, ölçüldü ≈%99.9), sonraki adımlar sönümler.
- **Pore-akışkan rijitliği Kw/n:** gerçek su bulk modülü Kw=2.0e6 kPa / porozite n=e/(1+e) (Verruijt; Sci.Man §4) →
  neredeyse sıkışmaz → cv ≈ k·Eoed/γw. v1 tek temsilci porozite (çok-malzemede yaklaşık).
- **Drenaj sınırı:** Flow conditions FlowBC'den türetilir (Head/Seepage kenar = açık/drenajlı p=0; Closed/tanımsız
  = geçirimsiz doğal). Hiç akış BC'si yoksa model üstü açık (varsayılan); pasif (kazılmış) elemanlara dokunan
  düğümler de drene (pore DOF taşımaz). `active` maskesi çekirdeğe geçer → kazı/dolgu sonrası konsolidasyon doğru.
- **Zaman:** `Phase.duration` [gün] + `Phase.time_steps` (eşit adım, dt = duration/steps).
- **Nihai durum:** efektif gerilme σ' = init + D·B·v_final (LE iskelet; `recover_consolidation_stress`); tam
  sönümde p→0 → toplam = efektif. committed sonraki faza taşınır (SumMstage zinciri). Çıktı: oturma-zaman +
  fazlalık-pore-zaman eğrisi (GUI PlotLines = PLAXIS U-t/Curves) + nihai alanlar.
- **Çözücü (TAM, sparse):** birleşik saddle-point matrisi A=[K L; Lᵀ −(ΔtH+S)] **simetrik-indefinite**
  → SparseMatrixBuilder'a TAM toplanır (PARDISO sarmalayıcı üst üçgeni kendisi çıkarır), **PARDISO mtype=−2
  ile BİR KEZ faktörlenir** (sabit dt) + her adımda yalnız geri-yerine-koyma (factor-once-solve-many,
  stateful shared_ptr çözücü; motorun geri kalanıyla tutarlı). H ayrıca seyrek tutulur (RHS ΔtH·pₙ çarpımı).
  Çekirdek MKL'den BAĞIMSIZ: `ConsolidationSolveFactory` callback (build_problem'de kurulur); boşsa dense
  Eigen LU (CSR→dense) = MKL'siz referans yol (test_consolidation). **Sparse=dense BİREBİR doğrulandı**
  (Terzaghi sayıları özdeş). 8000-bilinmeyen guard'ı KALDIRILDI (sparse ölçekler — PARDISO motor genelinde).
- **Doğrulama `test_consolidation_gui`:** 1D Terzaghi kolonu GUI yolundan (Project→mesh→K0 fazı→consolidation
  fazı sürşarjla); undrained pore üretimi ≈q (%0.1), Tv≈2'de pore→0.10 kPa, U_FE(Tv) Terzaghi'ye **<%1.2**,
  nihai oturma qH/Eoed'e **%0.6**.
- **v1 dışı (dürüst guard):** MC/HS elastoplastik konsolidasyon GUI (§8 çekirdek hazır, GUI wiring sıradaki),
  axisym, konsolidasyonda yapısal elemanlar, "min pore"/"degree of consolidation" yük tipleri.

## 8. ELASTOPLASTİK (MC/HS) konsolidasyon — kuplajlı Newton (TAM: çekirdek + GUI)
LE çekirdeğin (§4-6) doğrusal-olmayan genellemesi. Efektif gerilme bünye return-mapping'inden (integrate_point,
solve_nonlinear ile AYNI doğrulanmış fonksiyon) gelir; teğet K_T duruma bağlı. Her zaman adımında **MONOLİTİK
kuplajlı Newton:**
```
[K_T  L ] [δv]   [ r_u ]   r_u = Δf − Δf_int(Δv) − L·Δp        (denge artımı, Δf_int = f_int(trial)−f_int(committed_n))
[Lᵀ  −S*] [δp] = [ r_p ],  r_p = ΔtH·p_n − (Lᵀ Δv − S* Δp)     (süreklilik; S* = ΔtH + S)
```
yakınsayana dek iterasyon; commit'te trial Gauss durumları yeni committed olur (patika-bağımlı). **Pore = AÇIK DOF**
(undrained wrapper YOK); su sıkışabilirliği S = ∫(n/Kw)NᵀN'de → t=0+ near-incompressible ⇒ undrained plastik tepki
kendiliğinden. K_T non-associated plastisitede NONSİMETRİK → A tam nonsimetrik → çözücü **PARDISO RealNonsymmetric
(mtype=11)**, her iterasyonda yeniden faktörlenir (factory). Sıcak `solve_nonlinear` döngüsü DEĞİŞTİRİLMEDİ —
return-mapping izole `solve_consolidation_plastic`'te çağrıldı (perf bit-birebir korundu).
- **Doğrulama (`test_consolidation_plastic`):** (1) **LE-İNDİRGEME**: elastoplastik(LE iskelet) çözücü, doğrulanmış
  LE konsolidasyon çekirdeğine **1.188e-16 (makine hassasiyeti)** BİREBİR indirgenir + Terzaghi U-Tv'ye <%2.3 →
  kuplajlı-Newton makinesi (montaj, residual, kuplajlı çözüm) KESİN doğru. (2) **MC (akma yok) == LE 0.0 kesin** →
  MC return-mapping dalı kuplajlı montajda doğru bağlı. (3) **AKAN-VAKA DRENAJLI-LİMİT (gold oracle)**: MC footing
  (φ=20°, kayma plastisitesi devrede) undrained yüklenip t→∞'a konsolide edilince efektif oturma **doğrudan
  drenajlı elastoplastik solve_nonlinear'a %0.3** (her adım yakınsadı, fazlalık pore→0.04). İÇGÖRÜ: φ>0'da undrained
  yük boşluk suyunda → efektif gerilme≈değişmez → undrained plastisite MİNİMAL; plastisite konsolidasyon boyunca
  KADEMELİ → her adım modest artım → doğal iyi yakınsama + küçük patika-bağımlılık.
- **GUI (TAM):** build_problem consolidation fazı aktif malzeme MC/HS ise (`nonlinear_soil`) plastik yolu çağırır
  (RealNonsymmetric factory; önceki fazın committed efektif Gauss durumları = initial_state); LE-only → LE yolu.
  test_consolidation_gui MC oedometer (K0 confined = elastik) GUI'den = Terzaghi BİREBİR (plastik wiring + doğru
  indirgeme teyidi). KALAN: HS undrained-B su-cap konsolidasyon kalibrasyonu (niceliksel), axisym, "min pore".

İlgili: [[effective-stress-formulation]] (undrained Kw/n), [[seepage-formulation]] (H matrisi + FlowBC drenaj),
[[plaxis-gap-analysis]] Faz A.2, [[material-model-architecture]] (integrate_point), [[hardening-soil-formulation]].
