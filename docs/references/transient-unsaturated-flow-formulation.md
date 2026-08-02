# Transient + Doymamış Akış & Fully-Coupled Flow-Deformation — Formülasyon

"Su hikayesini bitir" izi. KATAI'de **steady-state seepage** (confined/unconfined/seepage-face) + **Biot
konsolidasyon** (zaman-bağımlı, doygun, sabit-k, elastoplastik kuplajlı-Newton) zaten TAM. Bu doküman
kalan üç boşluğu PLAXIS'in kendi denklemleriyle kilitler:

1. **Transient (zaman-bağımlı) yeraltısuyu akışı** — deformasyonsuz, depolama terimli.
2. **Doymamış akış** — van Genuchten SWCC + Mualem bağıl geçirgenlik (Picard altyapısı zaten var).
3. **Fully coupled flow-deformation** — af67414 konsolidasyon çekirdeğinin doygunluk + zamanla-değişen-BC
   genellemesi (PLAXIS'in en genel analizi).

**Kaynaklar (kilitli):**
- **PLAXIS 2D 2025.1 Scientific Manual §3** (Groundwater flow theory, Eq 3-1…3-37) — birebir okundu;
  §3.1.2 süreklilik (Eq 3-8/3-10), §3.3 FE ayrıklaştırma (Eq 3-31…3-37). **Fully-coupled** süreklilik
  = §3.1.2 Eq 3-8 (Sci.Man §4 konsolidasyonun doygunluk-terimli genellemesi).
- **PLAXIS 2D Reference Manual** — §7.4.4 fully-coupled flow-deformation analiz tipi; akış BC tipleri.
- **PLAXIS 2D Material Models Manual** — van Genuchten/Mualem (g_a, g_n, g_l) + approximate (Hypres/Staring)
  veri setleri. *(Exact parametre eşleşmesi: implementasyon oturumunun ilk işi — MMM PDF'i.)*
- **van Genuchten (1980)**, *A closed-form equation for predicting the hydraulic conductivity of
  unsaturated soils*, SSSAJ 44:892-898 — SWCC kapalı formu.
- **Mualem (1976)**, *A new model for predicting the hydraulic conductivity of unsaturated porous media*,
  WRR 12:513-522 — bağıl geçirgenlik.
- **Bishop (1959)** — doymamış efektif gerilme σ' = σ − χ·u_w·m (χ = S_eff).
- **Lewis & Schrefler**, *The Finite Element Method in the Static and Dynamic Deformation and
  Consolidation of Porous Media* (2. baskı) — doygun/doymamış kuplajlı FE referans kitabı.
- **Liakopoulos (1965)** — drenaj kolonu deneyi = doymamış kuplajlı altın-benchmark (çok-kod).
- **Carslaw & Jaeger** (ısı iletimi = lineer difüzyon) / **Theis (1935)** / **Bear (1972)** — transient
  akışın kapalı-form oracle'ları.

Konvansiyon: tension-pozitif (basınç p_w < 0 emme/suction), tüm kod tabanıyla ve [[consolidation-formulation]],
[[seepage-formulation]] ile aynı. p_w = boşluk suyu basıncı, head h = y + p_w/γ_w.

---

## 1. Transient akış (Sci.Man §3.1.1–3.1.2)

Darcy (Eq 3-7): q = (k_rel/γ_w) kᵃᵗ (∇p_w + ρ_w g),  k = k_rel·kˢᵃᵗ (Eq 3-4; k_rel = doygunluğa bağlı oran).

**Genel süreklilik (Eq 3-8 — fully-coupled, doygunluk-terimli):**
```
∇ᵀ(ρ_w q) + S·mᵀ ∂ε/∂t − n(S/K_w − ∂S/∂p_w) ∂p_w/∂t = 0        (3-8)
```
- `S·mᵀ ∂ε/∂t` = deformasyon kuplaj terimi (iskelet hacim değişimi). **Transient akış-only**'de katı
  deplasmanları ihmal → bu terim düşer (Eq 3-10).
- `n·S/K_w` = su sıkışabilirliği (doygun konsolidasyon S matrisi); `n·∂S/∂p_w` = **doymamış nem-kapasitesi**
  (SWCC eğiminden; doygun bölgede ∂S/∂p_w=0 ⇒ konsolidasyona indirgenir).

**Transient akış-only süreklilik (Eq 3-10):**
```
∇ᵀ(ρ_w q) − n(S/K_w − ∂S/∂p_w) ∂p_w/∂t = 0                      (3-10)
```
Steady-state (∂p_w/∂t=0) → ∇ᵀ(ρ_w q)=0 (Eq 3-11) = mevcut seepage çekirdeği.

### 1.1 FE ayrıklaştırma (Sci.Man §3.3, Eq 3-31…3-37)
p_w = N·p_wₙ (Eq 3-27, mevcut tri6/tri15 N). Galerkin → yarı-ayrık sistem (Eq 3-31):
```
−H·p_wₙ − S·(dp_wₙ/dt) = q_p
```
```
H   = ∫ (∇N)ᵀ (k_rel/γ_w) kˢᵃᵗ (∇N) dV                          (3-32)  [permeabilite — seepage He'nin k_rel·γ_w hali]
S   = ∫ Nᵀ ( nS/K_w − n ∂S/∂p_w ) N dV                          (3-33)  [depolama/sıkışabilirlik + nem-kapasitesi]
q_p = ∫ (∇N)ᵀ (k_rel/γ_w) kˢᵃᵗ ρ_w g dV  −  ∫ Nᵀ q̄ dΓ            (3-34)  [yerçekimi + verilen-akı]
```
**K_w/n (Eq 3-35):** = 3(ν_u−ν)/((1−2ν_u)(1+ν))·K_skeleton, ν_u default 0.495 — **konsolidasyon docu §2 ile
BİREBİR aynı** (yeni-açılan malzemede K_w ihmal).

**Zaman integrasyonu (Eq 3-36, α=1 tam-implicit; PLAXIS varsayılanı):**
```
−(αΔtH + S)·Δp_wₙ = ΔtH·p_wₙ⁰ + Δt·q_p          (α=1)          (3-36)
```
⇒ **bu, KATAI konsolidasyon saddle-point'inin `S* = αΔtH + S` sağ-alt bloğunun ta kendisi** (mekanik
kuplaj L bloğu çıkarılmış). `analysis/consolidation.hpp` H ve S'i zaten kuruyor + backward-Euler factor-once
yapıyor → **transient akış-only ≈ mevcut montajın mekaniksiz dalı.** Steady-state limiti Eq 3-37.

### 1.2 Sınır koşulları (Sci.Man §3.2)
Mevcut (seepage): Closed q_n=0 doğal (3-16) ✅ · Inflow/Outflow q̄ (3-17/18) ✅ · Head h=h̄ (3-19) ✅ ·
Seepage/phreatic (3-22) ✅. **Yeni (transient):**
- **Zamanla-değişen head** (Eq 3-23): h̄ = h̄(t); fazda zaman-tablosu. **rapid drawdown** için kritik.
- **Infiltration** (Eq 3-21): ponding/drying karışık BC (q̄ akı, kapasite aşılırsa h=z+h̄_p,max'a geçer,
  kuruyunca h=z+h̄_p,min). Aktif-set (seepage-face altyapısıyla aynı desen).
- **Wells** (Eq 3-24/25): nokta kaynak/kuyu Q=±Q̄ (iç düğüm source); **Drain** (Eq 3-26): iç seepage BC.

---

## 2. Doymamış akış — van Genuchten SWCC + Mualem k_rel

S(p_w) ve k_rel(p_w) doygunluk-su-karakteristik eğrisinden (SWCC). PLAXIS van Genuchten modelini kullanır
(parametreler g_a [1/m], g_n, g_l). Emme ψ = −p_w/γ_w (p_w<0 ⇒ ψ>0):

**Effective saturation (van Genuchten 1980):**
```
S_e = (S − S_res)/(S_sat − S_res) = [ 1 + (g_a·ψ)^{g_n} ]^{−g_m} ,   g_m = 1 − 1/g_n   (ψ>0; ψ≤0 ⇒ S_e=1)
```
**Bağıl geçirgenlik (Mualem 1976):**
```
k_rel = S_e^{g_l} · [ 1 − (1 − S_e^{1/g_m})^{g_m} ]^2          (g_l = Mualem pore-connectivity, default 0.5)
```
**Nem kapasitesi** (Eq 3-33 S matrisi için): ∂S/∂p_w = (dS/dψ)·(−1/γ_w), dS/dψ van Genuchten'in
analitik türevi (kapalı-form; lokal Newton/Picard için gerekli).

**Picard (Sci.Man §3.3): doymamış set yüksek-nonlineer** → her zaman adımında k_rel(p_w), S(p_w) güncellenir
→ lineer set implicit adımda çözülür → yinele. KATAI'de `solve_unconfined_seepage` zaten Picard + under-relax
(relax≈0.15) yapıyor; tek değişiklik kaba lineer-transition k_rel'i **van Genuchten/Mualem ile değiştirmek**.
NOT (Eq 3-15 dipnotu): hidrolik gradyan k_rel<0.99'da 0 alınır (yalnız doygun hacimde tanımlı).

---

## 3. Fully-coupled flow-deformation (Ref §7.4.4)

Konsolidasyon (Sci.Man §4, doygun) + doymamış + zamanla-değişen hidrolik BC. Monolitik sistem
**af67414'te ZATEN var**:
```
[K_T  L ] [δv]   [ r_u ]
[Lᵀ  −S*] [δp] = [ r_p ]                         S* = αΔtH + S
```
Genelleme (3 ekleme):
1. **Doygunluk kuplajı:** süreklilik Eq 3-8'deki `S·mᵀ∂ε/∂t` (doygun konsolidasyonda S=1) → L bloğu S ile
   ölçeklenir; S, ∂S/∂p_w doymamış bölgede SWCC'den (§2).
2. **Doymamış efektif gerilme (Bishop):** σ' = σ − χ·p_w·m, **χ = S_eff** (PLAXIS Bishop'ı). Doygun
   bölgede χ=1 ⇒ klasik Terzaghi (mevcut). Bu, `integrate_point`'in gördüğü efektif gerilmeyi etkiler.
3. **Zamanla-değişen BC:** fazda hidrolik BC (head/akı) zaman-tablosu → q_p her adımda güncellenir.
   (Mevcut konsolidasyon yükü yalnız t=0+ SumMstage; burada akış BC'si de zamanla değişir.)

Doygun + sabit-BC limiti **bit-birebir** mevcut konsolidasyona inmeli (regresyon).

---

## 4. Uygulama merdiveni (basitten karmaşığa, her adım izole doğrulanır)

| Adım | İçerik | Mevcut kod kaldıracı | Oracle (referans) | Hedef |
|---|---|---|---|---|
| **W1 ✅** | Transient akış-only (depolama+backward-Euler, zamanla-değişen head/akı BC) | `consolidation.hpp` H,S montajı + factor-once; seepage BC | 1D lineer difüzyon erfc / ani-head step (= Terzaghi matı, Carslaw-Jaeger/Bear); Theis transient kuyu | <%2 |
| **W2 ✅** | Doymamış: van Genuchten S_e(ψ)+Mualem k_rel + nem-kapasitesi | `solve_unconfined_seepage` Picard + under-relax | (a) SWCC nokta-değerleri formülle round-off; (b) Richards infiltration (Celia 1990 / Vauclin 1979) | <%5 |
| **W3 ✅** (LE) | Fully-coupled (doygunluk kuplajı + Bishop χ=S_eff + zaman-BC) | af67414 monolitik kuplajlı-Newton | (a) Terzaghi doygun-limit **bit-birebir** regresyon ✅ (max\|Δ\|=0); (b) **Liakopoulos** drenaj (çok-kod bant) [takip]; (c) rapid drawdown [takip] | bant-orta |
| **W4 ✅** | GUI: PhaseType::TransientFlow + FullyCoupled; zaman-BC; head/pore/oturma + doygunluk alanı | `PhaseType` deseni + PlotLines | GUI-yol = çekirdek (test_water_gui) | <%5 |
| **W3p ✅** (MC/HS) | Elastoplastik fully-coupled (W3 doymamış terimleri + return-mapping iskelet) | `consolidation_plastic` coupled-Newton + W3 katsayıları | (a) doygun-limit **bit-birebir** = consolidation_plastic ✅ (LE + MC-akan, max\|Δ\|=0); (b) LE-doymamış = W3 sabit-nokta (~%1e-5) ✅ | bit-/tol |

**Doğrulama disiplini:** matematik (yukarısı, PLAXIS §3 birebir) → izole test → kapalı-form (<%1-2) →
kuplajlı → çok-kod bant ortası. Kredibilite paralel: Terzaghi ✅ + **Liakopoulos** + **rapid drawdown** =
geoteknikteki en güçlü kanonik çok-kod su-benchmark'ları (band-match kanıtı için ideal).

---

## 5. Uygulama durumu

### W1 — TAM ✅ (`analysis/transient_flow.hpp`, `test_transient_flow`)
Head formunda transient doygun akış: `S·(dh/dt) + H·h = q`, backward-Euler (α=1), `A=(S+ΔtH)` bir kez
faktörlenir (factor-once-solve-many; A SPD → PARDISO mtype=2, boşsa dense LU referans). H = seepage
`element_conductivity`'nin ta kendisi; S = ∫NᵀS_sN (özgül depolama S_s=n·γ_w/K_w). Zamanla-değişen Dirichlet
head `HeadBoundary(node,t)` — sınır KÜMESİ sabit (yalnız değer değişir) → factor-once korunur, lift her adım
RHS'e taşınır. Eleman-generic (tri6/tri15). **Doğrulama (3 kapalı-form oracle, hepsi hedefin çok altında):**
| Oracle | KATAI hata | Referans |
|---|---|---|
| Terzaghi izokron u(Z,Tv)/u0 | ~%1 | Fourier serisi Σ(2/M)sin(MZ)exp(−M²Tv) |
| Ani-head adımı h(d,t)/hs | <%0.4 | Carslaw-Jaeger erfc(d/2√(Dt)) |
| Kararlı-hale yakınsama | 9e-14 (bit-yakın) | lineer head = steady seepage |

D = k/S_s (difüzyon katsayısı). Steady limiti `assemble_seepage`'e bit-yakın iner ⇒ regresyon güvencesi.
**GUI'ye bağlı DEĞİL** (W4'te tüm su hikayesiyle birlikte gelecek; PhaseType::TransientFlow).

### W2 — TAM ✅
**W2a (`materials/water_retention.hpp`, `test_water_retention`):** van Genuchten (1980) S_e + Mualem (1976)
k_rel + analitik nem-kapasitesi dS/dψ. Birim-doğrulandı: doygun limit, monotonluk, asimptotik eğim −(g_n−1)
(kesin), Mualem n=2 kapalı-form çapraz-kontrol, ve **analitik kapasite = merkezi sonlu-fark ~1e-11**
(makine hassasiyeti). PLAXIS g_a/g_n/g_l konvansiyonu (van Genuchten modeli = literatür standardı; PLAXIS de
aynısını uygular).

**W2b (`solve_transient_unsaturated_flow`, `test_unsaturated_flow`):** transient doymamış Richards akışı, head
formu, **kütle-korumalı modified-Picard (Celia 1990)** — mixed (θ-tabanlı) form: `[M_C/Δt + H^m]δ = q − H^m h^m
− (1/Δt)[∫N(θ^m−θ^n) + M_ss(h^m−h^n)]`. **CONSISTENT kütle** (Gauss noktasında θ/k_rel/kapasitans) — lumping
YOK (tri6 köşe ∫N_corner=0 dejenerasyonunu önler), global kütle-korunumu 1ᵀ∫N(θⁿ⁺¹−θⁿ)=∫Δθ ile korunur.
Doğrulama: doygun-limit Terzaghi izokron ~%1; **global kütle-korunumu oran=1.000000** (kesin, Celia kriteri);
**Philip √t sorptivite** I(8)/I(2)=2.20 (√t-yasası 2.0; lineer 4.0 olurdu). Picard keskin ıslanma-cephesinde
~30 iter (kuru zemine infiltrasyon, doğal stiffness).

### W3 — TAM ✅ (LE iskelet) (`analysis/coupled_flow_deformation.hpp`, `test_coupled_flow`)
Fully-coupled flow-deformation (PLAXIS'in en genel analizi), LE iskelet. Biot konsolidasyonun (consolidation.hpp)
doymamış genellemesi: saddle `[K  L_χ; L_χᵀ  −(ΔtH_kr+S_st)]`, **3 ekleme** (saturated limitte hepsi kaybolur):
(1) **Doygunluk kuplajı** L_χ=∫Bᵀ(S_eff·m)N (Eq 3-8 S·mᵀ∂ε/∂t); (2) **Bishop efektif gerilme** σ=σ'+χ·p·m,
χ=S_eff (LE'de χ yalnız L'yi etkiler, σ'=Dε); (3) **Doymamış depolama+geçirgenlik** S_st=∫Nᵀ(S·n/Kw+n·dS/dp_w)N,
H_kr=∫Gᵀ(k_rel·k/γw)G. emme ψ=−p/γw (p≥0⇒ψ≤0⇒doygun). Doymamış→nonlineer→adım başına PICARD (katsayıları lag).
DOĞRULAMA (test_coupled_flow): **(a) doygun-limit W3==solve_consolidation max|Δu|=max|Δp|=0 (KESİN bit-birebir)
→ tüm kuplaj makinesi doğru**; (b) Terzaghi U(Tv) ~%2; (c) **doymamış aktivasyon: başlangıç emme→S_eff=0.8575
(van Genuchten kapalı-form birebir), drenajla resaturasyon→1.0000, Picard 20 iter, S_eff∈(0,1]**.

### W4 — TAM ✅ (GUI entegrasyonu)
`PhaseType::TransientFlow` + `PhaseType::FullyCoupled` (model + project_io round-trip); build_problem branch'leri
(`InitialPhase::TransientFlow` → `solve_transient_unsaturated_flow`, flow-only; `FullyCoupled` →
`solve_coupled_flow_deformation`, LE iskelet); `flow_drained_nodes`/`flow_head_nodes` BC yardımcıları.
Malzeme `gw_ga/gw_gn/gw_gl/gw_Sres` (van Genuchten Groundwater sekmesi). GUI: faz tipi combo (5 tip), zaman
interval/steps, doygunluk alanı (Field::Saturation kontur), zaman-eğrileri (oturma + pore + min S), rapor
(txt+html) faz adları. **KRİTİK İŞARET DERSİ: kod tabanı TENSION-POZİTİF → sıkışma boşluk basıncı (yük altı) =
NEGATİF p → doygun; ψ_suction = +p/γw (NOT −p/γw); ilk sürümde ters işaret yükü her yerde van Genuchten'e sokup
Picard'ı yakınsamattı.** Doğrulama (test_water_gui): FullyCoupled GUI yolu = Terzaghi U(Tv) <%5 + doygunluk~1;
TransientFlow GUI yolu (basınçlandırma) pore gelişir + resaturasyon. GUI Picard tol 1e-7 (PARDISO re-factor için).

### W3p — TAM ✅ (elastoplastik MC/HS iskelet) (`analysis/coupled_flow_deformation.hpp` `solve_coupled_flow_deformation_plastic`, `test_coupled_flow_plastic`)
W3'ün (LE) ve elastoplastik konsolidasyonun (`consolidation_plastic`) **birleşimi** = PLAXIS 2D'nin EN GENEL
analizi: doymamış (van Genuchten/Mualem + Bishop χ=S_eff) Biot konsolidasyonu **MC/HS plastik iskeletle**.
İskelet σ' return-mapping'inden (`integrate_point`) gelir, teğet K_T duruma bağlı; akış katsayıları
(S_eff, k_rel, depolama) pore'a bağlı. Her zaman adımında **monolitik Newton-Picard**:
`[K_T  L_χ; L_χᵀ  −(ΔtH_kr+S_st)] [δv;δp] = r`, `r_u = Δf − (f_int(Δv)−Bbase) − L_χ·Δp`,
`r_p = ΔtH_kr·p_n − (L_χᵀΔv − (ΔtH_kr+S_st)Δp)`. K_T plastik return-mapping'den (Newton); L_χ/H_kr/S_st en son
pore tahmininde değerlendirilir (Picard-lag). K_T non-assoc'ta NONSİMETRİK → solve_factory **RealNonsymmetric
(mtype=11)**. Teğet modu consolidation_plastic ile aynı (HS→continuum, MC/LE→consistent). **DOĞRULAMA
(test_coupled_flow_plastic):** (1a) doygun-limit LE kolon = `solve_consolidation_plastic` **bit-birebir**
(max|Δu|=max|Δp|=0); (1b) doygun-limit **MC akan temel** (undrained→konsolide) = consolidation_plastic
**bit-birebir** (max|Δ|=0) → plastik return-mapping yolu kuplajla doğru bağlı; (2) **LE-doymamış** = W3
(`solve_coupled_flow_deformation`) AYNI sabit-nokta (~%1.6e-5; iki çözücü özdeş ayrık denklemleri çözer,
residual-tabanlı durdurma sert kw_over_n=1e9 sisteminde iterate-tabanlıyla makine hassasiyetine değil sub-%0.1'e
yakınsar); (3) **MC-doymamış** başlangıç emme drenajı → her adım yakınsar + admissible (S_eff∈(0,1]) +
resaturasyon, MC-akmasız = LE bit-birebir. **GUI'ye BAĞLI:** build_problem FullyCoupled dalı `nonlinear_soil`
(MC/HS) ise plastik çekirdeği çağırır (Consolidation deseni); LE yolu bit-birebir korunur. KRİTİK İŞARET DERSİ:
tension-pozitif → sıkışma boşluk basıncı NEGATİF p ⇒ doygun; suction ψ=+p/γw; doygun oracle u0=−10 kullanılmalı
(u0=+10 ψ=+1>0 yanlışlıkla suction → doymamış sanılır, regresyon patlar).

### KALAN — W3-takip sonrası
(a) **Liakopoulos** drenaj kolonu (çok-kod altın benchmark; yerçekimi-akışı gerektirir → head/total-pore formu,
dış deneysel veri); (b) **rapid drawdown** (PLAXIS/FLAC/GEO5 bandı); (c) FullyCoupled'da yapısal eleman desteği
(şu an soil-only); (d) test_water_gui'ye FullyCoupled+MC GUI-yol kantitatif kontrolü.

İlgili: [[consolidation-formulation]] (Biot çekirdek, S*=αΔtH+S, K_w/n), [[seepage-formulation]] (H matrisi,
Picard, BC), [[effective-stress-formulation]] (Bishop χ), [[plaxis-gap-analysis]] D/E/F, [[multi-code-validation-plan]].
