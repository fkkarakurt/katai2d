# Hardening Soil (HS) Modeli — Formülasyon (uygulama spesifikasyonu)

PLAXIS 2D'nin en çok kullanılan ileri bünye modeli. Gerilme-bağımlı rijitlik + iki akma
yüzeyi (kayma + cap) + MC kırılma. Felsefe: önce matematiği referansla kilitle, sonra en
küçük parçayı analitik doğrula (bkz `mohr-coulomb-formulation.md`).

**Kaynaklar (kilitli):**
- **Schanz, T., Vermeer, P.A., Bonnier, P.G. (1999).** *"The hardening soil model:
  formulation and verification."* Beyond 2000 in Computational Geotechnics (Plaxis
  Symposium, Amsterdam), Balkema, 281–296. **Birincil kaynak.**
- **Duncan, J.M. & Chang, C.-Y. (1970).** *"Nonlinear analysis of stress and strain in
  soils."* JSMFD ASCE 96(SM5):1629–1653. — hiperbolik gerilme-şekildeğiştirme temeli.
- **Rowe, P.W. (1962).** stress-dilatancy teorisi — mobilize dilatasyon.
- **Schanz & Vermeer (1996).** *"Angles of friction and dilatancy of sand."* Géotechnique
  46(1):145–151. — mobilize dilatasyon (cut-off, kritik durum).
- **PLAXIS Material Models Manual (Bentley)** — parametreler, p_ref, cap, hizalama referansı.
- (İleride) Benz (2007) — HS-small (küçük-şekildeğiştirme rijitliği).

## 0. Konvansiyon
Bu doküman **basınç-pozitif** (klasik triaxial/geoteknik): σ1 ≥ σ2 ≥ σ3 büyüklük, sıkışma
pozitif, q = σ1 − σ3 ≥ 0 deviatorik, p = (σ1+σ2+σ3)/3 ortalama. **Solver tension-pozitif**
olduğundan FE entegrasyonunda işaret dönüşümü gerekir (σ_solver = −σ_HS; bkz §7). Efektif
gerilme (HS daima efektif parametre c', φ' kullanır; su ayrı, bkz `effective-stress-formulation.md`).

## 1. Gerilme-bağımlı rijitlik (güç yasası)
Referans rijitlikler p_ref'te (genelde 100 kPa) verilir; gerçek rijitlik minör asal
efektif gerilmeyle ölçeklenir (m = üs: kumda ~0.5, kilde ~1):
```
E50  = E50_ref  · ( (c·cosφ + σ3·sinφ) / (c·cosφ + p_ref·sinφ) )^m     (sekant, %50 dayanım)
Eur  = Eur_ref  · ( (c·cosφ + σ3·sinφ) / (c·cosφ + p_ref·sinφ) )^m     (boşaltma-yeniden yükleme)
Eoed = Eoed_ref · ( (c·cosφ + σ1·sinφ) / (c·cosφ + p_ref·sinφ) )^m     (ödometre teğet; σ1'e bağlı)
```
Tipik: Eur_ref ≈ 3·E50_ref, Eoed_ref ≈ E50_ref (bağımsız belirlenir, sabit oran YOK).
ν_ur ≈ 0.2 (boşaltma-yeniden yükleme Poisson; elastik kısım). Eşdeğer biçim:
(c·cosφ + σ3·sinφ)/(...) = (c·cotφ + σ3)/(c·cotφ + p_ref) (PLAXIS'in c·cotφ biçimi).

## 2. Drained triaxial — hiperbolik primer yükleme (Duncan-Chang çekirdeği)
Sabit hücre basıncı σ3'te, deviatorik q artarken **primer** (ilk) yükleme eğrisi
hiperboliktir:
```
−ε1 = (1/Ei) · q / (1 − q/qa) ,      q < qf
```
- **Kırılma deviatoru (MC):** qf = 2(c·cosφ + σ3·sinφ)/(1 − sinφ)  [= σ3(Kp−1)+2c√Kp, Kp=(1+sinφ)/(1−sinφ)]
- **Asimptot:** qa = qf / Rf,   Rf = kırılma oranı ≈ 0.9 (qa > qf).
- **Başlangıç rijitliği:** Ei = 2·E50 / (2 − Rf).
- **%50 dayanımda sekant = E50 (özdeşlik):** q=qf/2'de sekant q/(−ε1) = Ei(1−Rf/2) = E50. *(Bu, E50 tanımının iç tutarlılık testidir — §6.)*
- **Kırılmada (q→qf):** mükemmel plastik (q sabit qf'te; deformasyon sürer). Eğri qf'te kesilir
  (hiperbol qa'ya gider ama fizik qf'te durur).
- **Boşaltma/yeniden yükleme:** elastik, dik eğim **Eur** (E50'den ~3× sert) → HS plastiktir,
  Duncan-Chang nonlineer-elastiğin aksine kalıcı şekildeğiştirme bırakır.

## 3. Kayma (shear/frictional) hardening akma yüzeyi
Sürtünme hardening'i deviatorik plastik şekildeğiştirmeyle (γ^p) genişler. Schanz (1999):
```
f_s = f̄(q, σ3) − γ^p ,      f̄ = (2/Ei)·q/(1 − q/qa) − 2q/Eur
```
- f̄ = toplam − elastik aksiyel şekildeğiştirme (hiperbol toplam − Eur elastik) ⇒ plastik kısım.
- **Hardening parametresi** γ^p ≈ −(2ε1^p − ε_v^p) (deviatorik plastik şekildeğiştirme; triaxial'de
  küçük dilatasyon için ≈ −2ε1^p).
- Yükleme yüzeyde (f_s=0, primer): γ^p = f̄(q) → §2 hiperbolünü üretir.
- **Teğet modül (primer triaxial):** Et = dq/d(−ε1) = Ei·(1 − q/qa)² (hiperbolün türevi).

## 4. Cap (volumetrik) hardening akma yüzeyi — OTORİTER FORM (Rocscience/PLAXIS HS doğrulama)
Sıkışmada (yüksek p) volumetrik plastik şekildeğiştirme; ön-konsolidasyon basıncı p_c hardening.
Eliptik cap (apeks q ekseninde):
```
F_c = (q*/α)² + p² − p_c² = 0                                              (Eq 15.11)
```
- **q\* özel gerilme ölçüsü (Lode-açısı bağımlı, Eq 15.12):**
  ```
  q* = q / f(θ) ,   f(θ) = (3 − sinφ) / (2(√3 cosθ − sinθ sinφ))
  ```
  Triaxial sıkışmada (σ2=σ3) f(θ) sabit → q* = q/sabit; bu sabit **α'ya emilir** (axisymmetric
  yollarda q̃=q yaklaşımı eşdeğer, K0'ı etkilemez). Genel gerilmede Lode bağımlılığı gerekir.
- **α = cap şekil faktörü.** DOĞRUDAN girdi DEĞİL — K0^NC'den türetilir (§4b).
- **Cap hardening GÜÇ YASASI (Eq 15.13, ÖNEMLİ — lineer DEĞİL):**
  ```
  ε_v^(p-cap) = (β/(1−m)) (p_c/p_ref)^(1−m)    ⇔    ṗ_c = (p_ref/β)(p_c/p_ref)^m
  ```
  β = cap hardening parametresi (DOĞRUDAN girdi DEĞİL — Eoed_ref'ten türetilir). Hardening
  modülü H_cap = dp_c/dε_v^(p-cap) = (p_ref/β)(p_c/p_ref)^m (gerilme-bağımlı). m=1: p_c=p_ref·exp(ε_v^pc/β).
  Hardening yalnız **cap'in ürettiği** volumetrik plastik şekildeğiştirmeyle.
- **Cap ASOSİYE** (akış yüzeye dik, dε^p=λ∂F_c/∂σ, dε_v^pc=λ·2p); kayma yüzeyi NON-ASOSİYE.

### 4a. p_limit (sayısal koruma)
σ3→0/negatif (kohezyonsuz çekme) stiffness denklemlerini (E50,Eur) sıfır/NaN yapar ve F_s
paydası patlar. **p_limit ≈ %10·p_ref**: minör asal gerilme bu limitin altına inmesin (σ3'ü
clamp et). PLAXIS/RS2-3 bunu yapar.

### 4b. Cap parametrelerinin (α, β) KALİBRASYONU — PLAXIS prosedürü (Eq 15.13 sonrası)
α ve β doğrudan girilmez; **ödometre testi simüle edilerek** K0^NC ve Eoed_ref'ten bulunur:
- Gerilme durumu K0 koşulunda, düşey gerilme = p_ref;
- Strain-kontrollü ödometre (aksiyel şekildeğiştirme artımı), **hem kayma hem cap aktif**;
- α, β öyle ki güncellenmiş gerilme durumu **K0^NC ve Eoed_ref üretsin** (2 koşul, 2 bilinmeyen).
PLAXIS bunu **sayısal** yapar (kapalı-form yok); biz de aynısını yaparız → gözlemlenen davranış
(K0^NC, Eoed) paritede.
- **BİLİNEN SINIRLAMA (literatür):** HS'in tahmin ettiği K0^NC, girilen K0^NC'den sapabilir —
  model bazı kullanıcılarca "kısıtlayıcı" bulunur (yüksek-K0 killer tam eşleşmeyebilir). **Bu
  PLAXIS'te de var** (model doğası, bug değil); doğru hedef PLAXIS'in *ulaştığı* davranış.

## 4c. DOĞRULAMA BENCHMARK'I — Berlin Sand III (PLAXIS 2014, Rocscience Tablo 15.1)
PLAXIS-doğrulanmış somut parametre seti (parite kanıtı):
```
p_ref=100 kPa, E50_ref=105 MPa, Eur_ref=315 MPa, Eoed_ref=105 MPa, m=0.55,
ν_ur=0.2, K0^NC=0.38, φ=38°, ψ=6°, c=1 kPa, Rf=0.9, T=0
```
Drained triaxial σ3=200: **qf=(c·cotφ+σ3)·2sinφ/(1−sinφ)=646 kPa** (Şekil 15.4 q-platosu ~640 ✓);
aksiyel gerilme platosu ~840 (=σ3+qf), volumetrik şekildeğiştirme dilatasyonu (ψ=6°). PLAXIS HS
ile RS2/3 HS bu eğrilerde örtüşür → hedefimiz aynı eğriler <%2-3.

## 4d. Deviatorik akma — OTORİTER (PLAXIS MMM Eq 6-8/6-9/6-10, TAM eşleşme)
PLAXIS MMM bizim formumuzu **birebir** kullanır (RS3'ün faktör-2'siz Eq 15.1 sadeleştirmesi değil):
```
f_s = f̄ − γ^p          (MMM Eq 6-8)
f̄ = (2/Ei) q/(1−q/qa) − 2q/Eur     (MMM Eq 6-9 — FAKTÖR-2 var, bizimle aynı)
γ^p = −(2ε1^p − ε_v^p) ≈ −2ε1^p    (MMM Eq 6-10 — hardening parametresi)
```
Bizim f̄ + γ^p (h_s=1−2R, §5) tam PLAXIS konvansiyonu. ψm=0'da hiperbolü üretir (test_hs_shear
rel<1e-9). qa, qf, Ei, Eur §1-2 ile aynı; q sınırı qf (MC kırılma).

## 4f. Cap deviatorik ölçüsü — von Mises q (SİMETRİK, kök-neden düzeltmesi)
Cap, **simetrik von Mises q** kullanmalı: `q² = 3J2 = ½[(σ1−σ2)²+(σ2−σ3)²+(σ3−σ1)²]`.
Eski δ·q̃ = σ1+(δ−1)σ2−δσ3 σ2,σ3'te ASİMETRİKti (∂/∂σ2=δ−1 ≠ ∂/∂σ3=−δ) → ödometrede σ2≠σ3
geliştirip K0'ı bozuyordu. 3J2 ve p² kuadratik ⇒ n_c=H_c·σ LİNEER kalır (closed-form return
korunur): **H_c = (3/α²)(I−⅓·11ᵀ) + (2/9)·11ᵀ**. (Rocscience q*=q/f(θ)'nin temeli; triaxial/
ödometrede f(θ) sabit α'ya emilir.) Bu düzeltme **kil K0=0.5931 vs 0.5933 TAM** yaptı (önce 0.43).

## 4g. Cap hardening anchor — KAPALI-FORM K_p (Itasca PH = PLAXIS-HS)
Cap plastik bulk modülü **K_p = K1·K2/(K1−K2)** (Itasca Plastic-Hardening dök., kamuya açık):
```
K1 = Eur_ref/(3(1−2ν))          (boşaltma bulk)
K2 = Eoed_ref(1+2K0nc)/3        (primer-yükleme bulk, K0nc dahil)
β = p_ref/(k·K_p),  k≈1 (ödometre düzeltmesi)
```
Cap hardening Eoed'e KAPALI-FORM bağlanır (`hs_cap_Kp`) → kalibrasyon: dış bisection α↔K0nc,
iç bisection k↔Eoed. Genel zemin (Eoed<E50) için **K0nc + Eoed İKİSİ DE ~TAM** (`test_hs_calibration`).

## 4e. K0^NC reprodüksiyonu — OTORİTER SONUÇ (PLAXIS MMM ile KAPANDI)
Berlin Sand III (Eoed_ref=E50_ref, çok sert cap, φ=38°) için girilen **K0^NC=0.38 ulaşılamaz**;
en-yakın-uygun ≈0.42. Bu bir solver kusuru DEĞİL — **PLAXIS'in belgelenmiş davranışıdır.** Bu
oturumda PLAXIS 2D Material Models Manual ve Brinkgreve (1994) kapalı-formları çekilip incelendi:

- **Cap yüzeyi (MMM Eq 6-26):** f_c = q̃²/M² + p'² − pp², q̃=σ1+(α−1)σ2−ασ3, α=(3+sinφ)/(3−sinφ).
- **M ↔ K0^NC (Brinkgreve 1994, MMM Eq 10-13 / Eq 6-29):** kapalı-form (yaklaşık M≈3.0−2.8·K0^NC).
- **Cap hardening ↔ Eoed (MMM Eq 6-27/6-28/6-30):** Ks/Kc ≈ (Eur_ref/Eoed_ref)·K0^NC/((1+2K0^NC)(1−2ν)),
  Ks_ref=Eur_ref/(3(1−2ν)); bizim güç-yasası β ile özdeş: β = p_ref·(Ks/Kc−1)/Ks_ref.
- **DOĞRUDAN ALINTI (MMM §6.4.3):** *"Depending on other parameters, such as E50, Eoed, Eur, there
  happens to be a certain range of valid K0^NC-values. K0^NC values outside this range are rejected
  by PLAXIS [...] the program shows the nearest possible value that will be used in the computations."*

**Erişilebilir küme kanıtı (sayısal, `study_hs_calibration`):** Berlin'de Eoed_ref=105000'i tutturmak
SERT cap ister → bu K0'ı ≥~0.42'ye iter; K0=0.38 ancak yumuşak capte (Eoed≈70000) elde edilir. (K0,Eoed)
2-serbestlik kuplajı (0.38, 105000)'i dışarıda bırakır — α/β taraması bunu doğruluyor. **q̃ (Eq 6-26)
asimetrik akışı ödometrede (de2=de3=0) σ2≠σ3 saçmalığı üretir; PLAXIS ödometre kalibrasyonu EKSENEL-
SİMETRİKtir, orada q̃ → von Mises'e indirgenir (σ2,σ3 ortalanır). Dolayısıyla bizim von Mises cap =
PLAXIS q̃ cap'in ödometre-yolundaki hâli; erişilebilir küme AYNI.**

**Sonuç:** kalibrasyonumuz **Eoed_ref'i TAM** tutturur (birincil rijitlik) ve **en-yakın-uygun K0'ı
(~0.42) raporlar** — tam olarak PLAXIS gibi. Genel zemin (Eoed<E50): K0^NC+Eoed İKİSİ DE ~TAM
(kil K0=0.5931 vs 0.5933). Berlin edge-case'i `test_hs_calibration` ile CI-korumalı (assert: Eoed tam,
K0 > 0.38 = restricted-range). Shear-baskın analizler (footing/şev/istinat) bu kısıttan ETKİLENMEZ.

## 5. Akış kuralı & mobilize dilatasyon (Rowe) — OTORİTER FORM (PLAXIS MMM §6, doğrulandı)
**Flow kuralı (PLAXIS MMM Eq 6-14):** `ε̇_v^p = sinψm · ε̇_q^p`. Tensör akış yönü (basınç-poz
asal, triaxial deviatorik ε_q^p=ε1^p−ε3^p): m_g=(1, R, R), R=ε3^p/ε1^p.
```
R = −(1 + sinψm)/(2 − sinψm)   ⇒   ε_v^p/ε_q^p = −sinψm   (DİLATASYON: basınç-poz ⇒ ε_v^p<0)
```
ψm=0'da R=−½ (volumetrik nötr). Hardening modülü (PLAXIS MMM **Eq 6-10**): γ^p=−(2ε1^p−ε_v^p) ⇒
dγ^p=h_s·dλ, **h_s = 1−2R = (4+sinψm)/(2−sinψm)** (ψm=0'da 2; hiperbol korunur). Eski sabit h_s=2
yalnız ψm=0'da doğruydu (Eq 6-10'un ε_v^p terimi atlanmıştı).

**KRİTİK İŞARET DÜZELTMESİ (2026-06-05, bu oturum):** önceki R=−(1−sinψm)/(2+sinψm) → ε_v^p/ε_q^p=
**+sinψm = KONTRAKSİYON** veriyordu (ψm>0). Dense kum (ψ=6°) kırılmaya yakın DİLATE etmeli (ε_v<0).
İşaret tersti; volumetrik hiç doğrulanmadığından (test_hs_berlin cap-OFF + volumetrik atlanmış)
yakalanmamıştı. Düzeltme PLAXIS Fig 15.4 volumetriğine karşı doğrulandı (`test_hs_berlin`).

**Kırılma platosu (perfectly-plastic MC, ÖNEMLİ):** q→qf'te hardening kayma yüzeyi f̄(q)−γ^p
qf'i geçemez (f̄ artar ama q clamp'lenir) → dormant kalır → dilatasyon durur, cap kontraksiyonu
baskın çıkar (YANLIŞ). PLAXIS'te kayma mekanizması kırılmada **mükemmel-plastik MC yüzeyine devreder**
(f=q−qf, yield grad (1,0,−1), hardening 0); akış (1,R,R) dilatant kalır → plato boyunca dilatasyon
sürer. `hs_integrate` bunu yapar (at_fail dalı) → Fig 15.4 net dilatasyonu (−0.005 @ %5) üretir.

Mobilize dilatasyon açısı ψ_m (Rowe stress-dilatancy, PLAXIS MMM Eq 6-15/6-16):
```
sinψ_m = (sinφ_m − sinφ_cs) / (1 − sinφ_m·sinφ_cs)          (φ_m mobilize sürtünme)
sinφ_cs = (sinφ − sinψ) / (1 − sinφ·sinψ)                   (kritik durum, φ_cs)
```
- Düşük mobilizasyonda (sinφ_m < ¾sinφ) ψ_m=0 (kontraksiyon/kritik-altı); üst sınır ψ_m≤ψ;
  dilatasyon cut-off (kritik boşluk oranı e_max'ta ψ_m→0, Eq 6-23/6-24 — şimdilik uygulanmadı).
- Volumetrik plastik şekildeğiştirme hızı: dε_v^p = sinψ_m · dγ^p.

## 6. Doğrulama planı (basitten karmaşığa)
1. **Stiffness güç yasası (P2.3a):** E50/Eur/Eoed gerilme bağımlılığı; σ3 oranında 2^m (c=0);
   E50 %50-sekant özdeşliği; qf = MC. *(Bu commit — `test_hardening_soil`.)*
2. **Drained triaxial hiperbol (P2.3a):** −ε1↔q eğrisi (§2), qf platosu, Eur boşaltma.
3. **Ödometre Eoed (P2.3b):** 1D sıkışma; gerilme-bağımlı Eoed, cap hardening; PLAXIS Eoed eğrisi.
4. **Tam HS tek-eleman (P2.3c):** triaxial+ödometre birlikte; Schanz (1999) doğrulama figürleri,
   PLAXIS HS örneğiyle <%5 (hatta daha sıkı — yüksek-derece eleman + tam entegrasyon).
5. **FE seviyesi (P2.3d):** multiaxial return mapping + tutarlı teğet; temel/şev/kazıda HS vs MC.

## 7. FE entegrasyonu — TAM İKİ-YÜZEYLİ (cap+shear, TAMAMLANDI)
- **`hs_integrate` (strain-driven, TAM iki-yüzeyli):** substepping + Koiter aktif-set (shear+cap),
  doğru Rowe dilatasyonu (§5) + kırılma platosu perfectly-plastic MC akışı. Kalibrasyon
  (`hs_oedometer_probe`/`hs_calibrate_cap`) ve tek-eleman triaxial (`test_hs_berlin`, Fig 15.4).
- **σ2=σ3 KORUNUMU (TUTARLI DRIFT, kritik):** shear drift düzeltmesi artık çıplak yield-gradyanı
  (1,0,−1) yerine **TUTARLI elastoplastik yön De·n_s boyunca** çeker (Potts & Gens 1985):
  δσ=−(f/(aᵀDe·n_s+h_s))·De·n_s. De·n_s (n_s=(1,R,R)) σ2,σ3'te simetrik ⇒ triaxial/oedometre
  kenarında (σ2=σ3, ör. axisym hoop=radial) σ2=σ3 KORUNUR. Çıplak gradyan yalnız σ1,σ3'ü itip σ2'yi
  bırakıyordu → her adım küçük asimetri birikip axisym oedometre K0'ını bozuyordu (0.53 vs 0.47).
  Düzeltme sonrası axisym oedometre K0=probe=kalibrasyon (`test_hs_axisym`), tüm tek-eleman testleri
  korundu (Berlin/calibration/integrate).
- **FE yolu (`integrate_point`→`hs_forward`→`hs_integrate`):** ROBUST. FE forward artık fragile
  stress-driven nested-Newton (`hs_return_principal`) yerine **strain-driven explicit-substepping
  integratörü `hs_integrate`** kullanır (kalibrasyon/oedometre-probe ile AYNI sağlam yol). hs_forward:
  full-Voigt elastik öngörücü (frozen Eur) → trial principal σ_tr (comp-pos desc) + committed principal
  σ_n; principal elastik şek.değ. artımı `deps_p = C_e·(σ_tr−σ_n)` (aynı frozen Eur ⇒ substep-1
  elastik trial'ı birebir kurtarır, elastik adım byte-identical) → `hs_integrate(σ_n, γ^p, pp, deps_p)` →
  coaxial reconstruction. **s3_stiff floor = 0.1·p_ref** (Voigt öngörücü = compliance = hs_integrate'in
  Eur'u tutarlı). σ_HS = −σ_solver.
- **ANALİTİK tutarlı teğet (FD DEĞİL):** eski sayısal merkezi-fark teğeti (6 ekstra entegrasyon/Gauss)
  substep aktif-set anahtarlamasından GÜRÜLTÜLÜ idi → global yakınsama rel~1e−3'te tabanlanıp footing
  BVP'sinde ıraksıyordu. Yeni: `hs_integrate`'in principal continuum teğeti D_pp'den ANALİTİK kurulur.
  Stress-to-stress principal Jacobian `J_comp = D_pp·C_e`; tension-poz sıralı konvansiyonda
  `J(i,j)=J_comp(2−i,2−j)` (işaret çevrimi Jacobian'ı değiştirmez, sıralama tersine döner). Bu J,
  **MC ile PAYLAŞILAN** `principal_consistent_tangent` (mohr_coulomb.hpp) yardımcısına verilir — eigenframe
  dönmesi (spin) + σ_zz kuplajı tek yerde (Sysala 2016 §5). MC closed-form region J ile, HS substepping
  continuum D_pp ile aynı makineyi kullanır. mc_return_mapping bu yardımcıya refactor edildi (byte-identical).
- **AXISYMMETRY TAM (`integrate_point_axisym`):** PLAXIS gibi HS hem plane strain hem axisymmetric. Çekirdek
  KİNEMATİK-AGNOSTİK (`hs_return_core`): (r,z) bloğu in-plane Mohr çemberi, hoop σ_θ üçüncü principal —
  plane strain'le AYNI yapı, principal return birebir aynı. Yalnız elastik öngörücü/operatör farklı: trial
  AXİSİMETRİK De(Eur) ile (hoop GERÇEK şek.değ.), tutarlı 4×4 teğet D_T=algo_jacobian·De_axisym (MC ile
  birebir aynı kalıp). `test_hs_axisym`: oedometre K0=probe=kalibrasyon + σ_r=σ_θ + dairesel temel BVP yakınsar.
- **Başlangıç durumu (FE init, ÖNEMLİ):** jeostatik K0 ön-gerilmesi admissible olmalı:
  `hs_initial_pp` (cap: pp=p_eq(σ0)·OCR, NC=cap üzerinde) + `hs_initial_gamma_p` (kayma: γ^p=f̄(q0),
  K0 durumu kayma yüzeyinde). Aksi hâlde ilk return büyük tutarsız düzeltme yapar → solver ıraksar.
- **BİRLEŞİK K0 PROSEDÜRÜ ŞART (plastik soil):** jeostatik denge `constant_force=∫BᵀΣ0` ile sabit tutulur,
  YALNIZ dış yükler ramp'lanır (gravity'yi ramp'lamak ara adımlarda sahte dengesizlik → düşük-confinement
  yüzey zemininde sahte plastisite → ıraksama). `assemble_internal_force(init)` kullan.
- **Yakınsama karakteri:** continuum teğet → LİNEER (kuadratik değil) global yakınsama; PLAXIS-gerçekçi
  tolerated-error (~%1=1e−2) + küçük yük artımları gerekir. build_problem (GUI) HS için steps=40, tol=1e−2
  (MC closed-form CONSISTENT teğet → quadratic, 1e−6). Continuum teğetin gerçek algoritmik teğete
  yükseltilmesi (kuadratik) gelecek incelik; near-cohesionless footing kenar tekilliği (c≈0, σ3→0) hâlâ
  chatter yapar (PLAXIS dahil bilinen FE patolojisi, küçük c veya mesh sıklaştırma ile çözülür).
- **Doğrulama:** (1) `test_hs_fe_oedometer` konfine HS kolon (birleşik K0 + surcharge): yakınsar (61/55 iter,
  eski fragile yol 169/175), cap AKTİF, cap-ON > cap-OFF oturma, K0≈0.39 (NC, Jaky 0.384). (2)
  **`test_hs_footing` (YENİ, kayma-baskın BVP):** şerit temel servis yükü (q=150, ~%38 kapasitesi);
  eski fragile yolda HANG/ıraksama, artık load_factor=1.0 + fiziksel oturma −2.96cm. Tek-eleman paritesi
  `test_hs_berlin`/`test_hs_calibration`.
- **UNDRAINED triaxial (PLAXIS Fig 15.5, `test_hs_undrained_triaxial`):** Berlin Sand III drenajsız triaxial,
  sabit-hacim (εv=0) efektif-gerilme yolu (u=σ3_hücre−σ3'). EŞLEŞEN: boşluk basıncı piki +92 vs PLAXIS ~+100;
  kontraksiyon→dilatasyon işaret dönüşü; q drenajlı qf'i (645) aşar (dilatant güçlenme, q@3%=878). **BİLİNEN GAP
  (cap-Lode kuplajı):** dilatasyon MAGNİTÜDÜ — cap-ON u_min=−84 vs PLAXIS −200 (cap-OFF −352 ile bracketlenir);
  **von Mises cap (q²=3J2) triaxial'de fazla kontraksiyon yapar.** PLAXIS Lode-bağımlı q̃=q/f(θ) (Eq 15.12, triaxial
  comp f(θ)≈0.66) kullanır → cap farklı angaje olur. von Mises cap önceden ödometre σ2=σ3'ünü bozduğu için seçilmişti;
  σ2=σ3 artık TUTARLI drift ile korunduğundan tam q̃-cap'e geçiş yolu açık (derin/riskli — yeniden kalibrasyon +
  drenajlı re-doğrulama).
- **UNDRAINED NC KİL (ψ=0, `test_hs_undrained_clay`) — TAM (von Mises cap EXACT):** ψ=0 ⇒ kayma dilatasyonu yok
  ⇒ cap-Lode inceltmesi GEREKSİZ (gap yalnız dilatant-undrained'de). NC kil drenajsız (sabit-hacim) yalnız
  KONTRAKSİYON yapar (u monoton +, p' düşer, ESP sola), efektif MC kırılma zarfında biter:
  q_f=M·p'+(6c cosφ)/(3−sinφ), M=6sinφ/(3−sinφ) → **rel hata %0.4** (PLAXIS Undrained-A NC-kil yolu, kesin).
- **UNDRAINED BVP (`test_hs_undrained_bvp`) — wrapper+HS:** global Undrained-A (D_u=D'+Kw/n·mmᵀ, total=σ'+Kw/n·εv·m)
  HS ile çalışır: şerit temel kil, εvol≈8e-6 (sıkışmaz ✓), boşluk basıncı gelişir ✓, undrained drenajlıdan ~4× rijit ✓.
  **KISIT:** yüksek undrained yük yakınsama TAVANI (~q_eff 9.4 kPa, geostatik ref + continuum-teğet tabanı, νu'dan
  BAĞIMSIZ) → continuum→ALGORİTMİK consistent teğet gerekir (drenajlı footing yavaşlığı + undrained tavanı AYNI kök).
  Modest yükte (q=6) yakınsar (lf=1, 53 iter). 73/73 test.

İlgili: [[literature-review]], efektif gerilme [[effective-stress-formulation]], MC
[[mohr-coulomb-formulation]].

## 9. Global teğet: HIBRIT continuum + SAYISAL CONSISTENT (Perez-Foguet & Rodriguez-Ferran & Huerta) - TAM (2026-06-13)

**Problem (linchpin):** hs_integrate'in continuum teğeti (son-aktif-set D_ep) global Newton'u LINEER yakinsatir;
duz/konfine rejimde adim basina ~3-4 iterasyonla yeterli AMA kayma-baskin dusuk-cevre-basinci rejiminde (temel
kenari, serbest yuzey) artimlar hic yakinsamaz: GUI yolu HS footing (unstructured mesh, K0+SumMstage)
load_factor=0.976 (coarse) / 0.844 (fine) tavaninda kalir.

**Cozum 1 - sayisal consistent teğet:** tam substep'li guncellemenin (aktif-set degisimleri + kirilma platosu +
drift duzeltme DAHIL) toplam-artima gore ileri-fark turevi: D(:,j) = [sigma(deps + h e_j) - sigma(deps)] / h,
h = sqrt(eps_makine) (1+|deps_j|). Kaynak: A. Perez-Foguet, A. Rodriguez-Ferran, A. Huerta,
"Numerical differentiation for local and global tangent operators in computational plasticity",
CMAME 189 (2000) 277-296; ve "Consistent tangent matrices for substepping schemes", CMAME 190 (2001) 4627-4647.
**KRITIK UYGULAMA DETAYI (2001 makalesinin ana pratigi):** perturbe kosular taban kosunun ALT-ADIM SAYISINA
sabitlenir (hs_integrate nsub_fixed) - serbest birakilirsa nsub=ceil(...) perturbasyonla n -> n+1 sicrar ve iki
kosu arasindaki entegrasyon-hatasi farki (~1e-2 kPa) FD sinyalinin (~E*h ~ 1e-3 kPa) ustune cikar, kolon bozulur.
Elastik adimda harita lineer (dondurulmus Eur) => continuum = consistent, FD kosulmaz. Plane strain 3 kolon,
axisym 4 kolon (perturbe predictor ayni dondurulmus De_ur ile kurulur, hs_return_core cagrisi nsub-sabit).

**Olcum (A/B, ayni build):** consistent teğet GUI tavanini COZER (0.976/0.844 -> 1.000/1.000) ama duz rejimde
~2.5x pahali (teğet montaji 4-5x: 3-4 ekstra hs_forward/Gauss) ve oradaki yakinsama zaten tol=1e-2'ye 3-4
iterasyonda ulasiyor (kuadratik fayda gorunmuyor; aktif-set/spektral-siralama kiviklari FD kolonlarini gecis
noktalarinda harmanlar - kNone/kContinuum/kConsistent yorumu, material_model.hpp).

**Cozum 2 - HIBRIT POLITIKA (nonlinear_solver):** artim CONTINUUM teğetle baslar; yakinsamazsa dlam yarilanmadan
once AYNI artim CONSISTENT teğetle yeniden denenir; artim commit olunca continuum'a geri donulur (artim-basina
yukselme). Yalniz HS malzeme varken devrede => LE/MC yolu BIREBIR eski davranis. Olcum (study_perf GUI yolu):
- HS footing coarse: 6.76s lf=0.976 (FAIL) -> 10.06s lf=1.000 (hibrit)  [saf consistent: 14.51s]
- HS footing fine:  14.87s lf=0.844 (FAIL) -> 16.06s lf=1.000 (hibrit)  [saf consistent: 41.77s]
Yani neredeyse continuum maliyetine TAM yakinsama; undrained yuksek-yuk tavani da ayni kokten cozulmeye aday
(bkz test_hs_undrained_bvp kisidi). Kalici-yukselme varyanti olculup REDDEDILDI (gecici takilma sonsuza dek
pahali modda birakiyor: structured footing 829 vs 594 iter).
