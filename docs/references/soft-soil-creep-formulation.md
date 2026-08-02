# Soft Soil Creep (SSC) Modeli — Kilitli Formülasyon (PLAXIS MMM §11, birebir)

**Birincil kaynak:** PLAXIS 2D 2025.1 Material Models Manual, Bölüm 11 "Soft Soil Creep model
(time dependent behaviour) [ADV]" (metin scratchpad çıkarımından; denklem numaraları manuale).
Kullanıcı ilkesi: bu modelin "doğrusu" tanım gereği PLAXIS'tir. Durum: **AŞAMA 1 (malzeme-noktası
çekirdeği) + AŞAMA 2 (FE entegrasyonu, §8) tamam** (`materials/soft_soil_creep.hpp` +
`material_model.hpp ssc_forward` + `test_soft_soil_creep`); GUI/proje-dosyası/faz-zamanı arayüzü +
Konsolidasyon-fazı V&V (MMM §17.4/17.5) = Aşama 3. İşaret kabulü çekirdekte basınç-POZİTİF (SS/HS
gibi; manual basınç-negatif yazar — çevrildi).

## 1. 1B creep zinciri (§11.2-11.4)
- Buisman (1936): sabit efektif gerilmede ε_c ∝ log(t) — ikincil sıkışma. Garlanger (1972) /
  Bjerrum (1967): e-log(t) biçimi; Butterfield (1979): logaritmik (Hencky) şekildeğiştirme.
- Janbu (1969) yöntemi: 1/ε̇ – t doğrusu → τ_c kesişimi ve c eğimi deneysel olarak belirlenir.
- **Temel kabul (Eq 11-11/12):** tüm inelastik şekildeğiştirme zamana bağlıdır (ani plastik bileşen
  YOK — kırılma hariç, §11.7) ve ön-konsolidasyon TAMAMEN birikmiş creep'ten büyür:
  **σ_p = σ_p0·exp(ε_c/(b−a))** (basınç-pozitif; b−a = 1B'de λ−κ karşılığı).
- **Diferansiyel yasa (Eq 11-19):** ε̇ = a·σ̇'/σ' + (c/τ)·(σ'/σ_p)^((b−a)/c). τ = **1 gün** —
  standart ödometrenin 24 saatlik kademe süresi NC çizgisinin TANIMIDIR (Eq 11-13/14: kademe sonu
  τ'de durum tam NC çizgisine oturur; τ_c ve t_c hakkında ek varsayım gerekmez).

## 2. 3B genişleme (§11.5)
- **Eşdeğer basınç (Eq 11-20):** p_eq, p–q̃ düzleminde Modified Cam-Clay (Roscoe & Burland 1968)
  elipsi üzerinde sabittir; q̃ HS/SS ile AYNI deviatorik ölçü (q̃ = σ1+(δ−1)σ2−δσ3):
  **p_eq = p′ + q̃²/(M²(p′ + c·cotφ))** — SS cap fonksiyonu f̄ ile ÖZDEŞ ölçü (kod tek kaynak:
  `ss_initial_pp`/cap_normal makinesi yeniden kullanılır).
- **Hacimsel creep hızı (Eq 11-23):**
  **ε̇_v^c = (μ*/τ)·(p_eq/p_p^eq)^β,  β = (λ*−κ*)/μ***
  p_p^eq = p_p0^eq·exp(ε_v^c/(λ*−κ*)) (üstel yaşlanma/pekleşme; Bjerrum "ageing").
- **Akış kuralı (Eq 11-26/30):** creep = zamana bağlı plastik şekildeğiştirme; potansiyel g = p_eq:
  **ε̇^c = (ε̇_v^c / (∂p_eq/∂p′)) · ∂p_eq/∂σ′** — hacimsel iz her zaman 1B'den türetilen hıza eşit
  (normalizasyon ∂p_eq/∂p′ = 1 − q̃²/(M²p̄²) ile). Elipsin kuru yanında (q̃ > M·p̄) bu terim
  negatiftir → creep hacimsel GENLEŞME/yumuşama üretir (modelin bilinen davranışı; MC §11.7 sınırlar).
- **M (Eq 11-37):** kullanıcı M girmez; K0NC'den Brinkgreve (1994) bağıntısıyla türetilir — SS ile
  AYNI formül (`M_from_K0nc`). (Eq 11-21'in M = f(φ_cv) yorumu bilgilendiricidir; kalibrasyon K0NC.)

## 3. Elastisite (§11.6) ve kırılma (§11.7)
- Elastik kısım SS ile aynıdır: K_ur = p′/κ* (Eq 11-32/33: E_ur girdisi YOK, κ*+ν_ur'den), üstel
  ln-yasa; ν_ur gerçek elastik sabittir (varsayılan 0.15; 1B boşaltmada K0 artışı ν/(1−ν) yönetir —
  κ*'ın 1B a/Cs ile KESİN bağıntısı yoktur, izotrop K0=1 yaklaşımı ile κ* ≈ 2a; §11.6 not).
- **Kırılma:** creep formülasyonu kırılma içermez → standart Mohr-Coulomb (c, φ, ψ) mükemmel-plastik
  olarak EKLENİR. Sıra (manual, §11.7): her gerilme noktasında ÖNCE creep güncellemesi, SONRA MC
  kontrol/düzeltme. Çekirdek bu sırayı birebir uygular (MC dönüşü tek kaynak: `softsoil::ss_step`
  pp=∞ ile — cap kapalı, aynı elastik yasa + kenar-kademeli MC).

## 4. Parametreler (§11.8) ve pratik notlar (§11.11)
λ* · κ* · **μ*** (modifiye creep indisi) · ν_ur (0.15) · c′ · φ′ · ψ · σ_t · K0NC (→M) ·
alternatif Cc/Cs/**Cα**/e_init: λ* = Cc/(2.3(1+e)), κ* ≈ 2Cs/(2.3(1+e)), **μ* = Cα/(2.3(1+e))**
(Tablo 11-3). Kaba kestirimler: μ* ≈ I_p(%)/500? — hayır: **λ* ≈ I_p(%)/500**, λ*/μ* = 15–25,
λ*/κ* = 2.5–7 (§11.8.1). Başlangıç p_p: OCR/POP'tan (§2.8; `ss_initial_pp` aynı p_eq ölçüsü).
**§11.11 uyarısı (ürünleşecek):** OCR=1.0 başlangıç gerçek dışı büyük başlangıç creep hızı verir
(hız ∝ OCR^−β, β~15-25); pratik öneri OCR ≥ 1.2–1.4. Dinamikte (§11.10) creep konturu içi çevrimler
yalnız elastiktir (histeretik sönüm yok → Rayleigh kullanılır) — SS ile aynı sınıf.

## 5. Kilit türetimler (çekirdek + V&V oracle'ları; buradan bağımsız doğrulanır)
1. **Sabit gerilme, NC (Buisman/Garlanger oracle):** p_eq sabit, p_p0 = p_eq:
   (p_eq/p_p)^β = exp(−β·ε_v^c/(λ*−κ*)) = exp(−ε_v^c/μ*) →
   dε/dt = (μ*/τ)e^(−ε/μ*) → **ε_v^c(t) = μ*·ln(1 + t/τ) KESİN** (ikincil eğim μ*/ln-çevrim;
   t≫τ'de Buisman). OCR'li başlangıçta: ε_v^c(t) = μ*·ln(OCR^−β·(t/τ)·β·μ*/(λ*−κ*)·...) yerine
   pratik pin: **başlangıç hızı oranı = OCR^−β** (aşağıda (b)).
2. **İzotropik gerilme RELAKSASYONU (deps = 0), NC:** dε_v = 0 → κ*·dp/p = −dε_v^c;
   p_p = p0·(p0/p)^(κ*/(λ*−κ*)) → p/p_p = (p/p0)^(λ*/(λ*−κ*)) → y = p/p0 için
   dy·y^(−(1+γ)) = −(μ*/(κ*τ))dt, γ = λ*/μ* → **p(t) = p0·(1 + (λ*/(κ*τ))·t)^(−μ*/λ*)** (kapalı
   form; V&V doğrudan pinler).
3. **τ'nin NC-çizgisi anlamı:** kademeli ödometre (her kademe 1 gün beklemeli) kademe sonlarında
   λ* NC çizgisi üzerinde ilerler (Eq 11-13/14) — V&V (f) bunu ölçer.
4. **Ölçü birliği:** p_eq ≡ SS f̄ → `ss_initial_pp` SSC için de başlangıç p_p verir; NC tohum
   p_eq = p_p (creep hızı tam (μ*/τ)).

## 6. Aşama-1 çekirdek şeması (`soft_soil_creep.hpp`) ve dürüstlük notları
- Şekildeğiştirme-sürümlü asal çekirdek + **zaman adımı dt**: (1) SS ile aynı üstel-kesin elastik
  öngörücü (p_tr = p_c·exp(Δε_v/κ*), G(p_tr)); (2) **creep düzeltmesi**: skaler L = Δε_v^c için
  geri-Euler: R(L) = L − dt·(μ*/τ)·(p_eq(σ(L))/p_p(L))^β = 0 — R monoton artan (L↑ ⇒ p_eq↓, p_p↑)
  → köşeli sekant+ikiye-bölme GARANTİLİ yakınsar; σ(L) üstel-ortalama dönüşle (iz-normalize m,
  trace(m)=1 ⇒ p(L) = p_tr·e^(−L/κ*)), yön trial'da DONUK (açık yön/kapalı büyüklük);
  (3) MC kontrolü: `softsoil::ss_step` pp=∞ (cap kapalı) — tek-kaynak MC + kenar kademesi;
  (4) alt-adım: SS gerinim ölçütü + creep büyüklüğü ölçütü (L_explicit/(0.05κ*)), nsub_fixed
  FD-teğet sabitlemesi için dışa raporlanır.
- DÜRÜST sınırlar: yön donukluğu ve elipsin kuru yanındaki ∂p_eq/∂p′ < 0 bölgesi (genleşmeli creep)
  alt-adımla sınırlanır; |∂p_eq/∂p′| tabanı (nose çevresi) sayısal korumadır; deviatorik creep yönü
  Aşama-1'de kapalı-formla PİNLENMEMİŞTİR (MMM §17.4 CRS/undrained-creep karşılaştırmaları FE
  aşamasının V&V'sine kayıtlı iş).

## 7. V&V (test_soft_soil_creep; oracle'lar §5 türetimleri — ölçülen sonuçlar)
(a) sabit izotropik gerilme: μ*·ln(1+t/τ) −%0.33 (t=2τ) / −%0.19 (t=100τ);
(b) OCR başlangıç hız oranı = OCR^−β: 65503 vs 2¹⁶ (−%0.05);
(c) izotropik relaksasyon kapalı formu +%0.01, yaşlanan p_p kesin;
(d) 24 saatlik kademeli ödometre: KARARLI kademe eğimi λ* +%0.2 (kademe-1 oturma geçişi hariç);
    kademe-sonu K0, K0NC'nin ÜSTÜNDE gezer (hızlı yüklemede elastik K0 düşer, beklemede creep
    toparlar — üst dönüm noktası) → M(K0NC) pini SÜREKLİ (CRS) yolda: K0 = 0.583 vs 0.577 (+%0.9;
    izotah çerçevesi oranı hızdan bağımsız kılar); CRS'te yanallar 1e-6 bağlı (köşe kademesi);
(e) hızlı drenajlı triaxial: MC çizgisi +%0.00.

## 8. Aşama 2 — FE entegrasyonu: zaman integrate_point'e nasıl girer
- `integrate_point(..., TangentMode, double dt_day = 0.0)`: kuyruk parametresi; 0 = creep yok
  (tüm eski çağıranlar bit-birebir; PLAXIS'te de zaman-aralıksız fazda SSC yalnız elastik+MC).
- STATİK çözücü: `NewtonOptions.time_interval` [gün] fazın süresi; her artıma **Δλ orantısıyla**
  paylaştırılır (dt = T·Δλ — PLAXIS'te zamanın SumMstage ile ilerlemesi; artım yarılanırsa zaman
  payı da yarılanır → fiziksel tutarlılık kendiliğinden). Montajcıda `InternalForceAssembler
  ::dt_day` üyesi; assemble adaptörü her artımda yazar.
- KONSOLİDASYON + tam-kuplajlı akış: kendi zaman adımları (gün) doğrudan geçer → creep +
  boşluk-suyu sönümü etkileşimi bünye düzeyinde kendiliğinden (faz-düzeyi V&V Aşama 3'te,
  MMM §17.5 deseninde). DİNAMİK: dt_saniye/86400 geçirilir (deprem süresinde creep ihmal
  sınıfı ama zaman muhasebesi sessizce sıfırlanmaz).
- SPEKTRAL İSKELET TEK KAYNAK: `ss_return_core` genelleştirildi (ssc işaretçisi + dt) — SS'in
  Voigt üstel-trial + rank-eşleme + koaksiyel rekonstrüksiyon kodu SSC için de birebir aynı
  (elastik yasa özdeş olduğundan tersinme geçerli); köşe/ulp sınıfı hatalar iki yerde yaşanmaz.
  FD teğet SABİT dt'de ve taban nsub'a kilitli (∂σ/∂ε|_t — Newton'un aradığı türev).
- **Aşama-2 V&V (g1-g3):** düzlem + eksenel integrate_point zinciri (şekildeğiştirme+zaman) ==
  çekirdek zinciri 1e-8; tri15 kolon BVP'si: zamanlı sıkıştırma + TUTULAN deplasmanla relaksasyon,
  her Gauss çekirdek zincirine ~2e-12, kolon 5 günde 223→146 kPa gerçek relaksasyon (zaman
  çözücüden uçtan uca akıyor).
