# HSsmall — Hardening Soil with Small-Strain Stiffness (Faz A.5)

Hardening Soil'in üstüne **küçük-şekildeğiştirme rijitliği**: çok küçük şekildeğiştirmelerde zemin çok
daha rijittir (G₀), artan kayma ile rijitlik degrade olur (Hardin-Drnevich). İnsanların PLAXIS'te ilk
sorduğu malzeme modeli. **Kaynak (kilitli):** PLAXIS 2D 2025.1 Material Models Manual §7 (Eq 7-1…7-12;
Benz 2006 small-strain overlay, Santos&Correia 2001, Hardin&Drnevich 1972, Masing 1926).

## 1. İki ek parametre (HS'in üstüne)
- **G₀^ref**: çok-küçük-şekildeğiştirme kayma modülü (ε<10⁻⁶, reference −σ'₃=p^ref). G₀^ref=E₀^ref/(2(1+ν_ur)) (Eq 7-12).
- **γ₀.₇**: sekant kayma modülü G_s'in 0.722·G₀'a düştüğü eşik kayma şekildeğiştirmesi (virgin loading).
Diğer tüm parametreler HS ile aynı.

## 2. Hardin-Drnevich (modifiye, Santos&Correia) — degradasyon yasası
```
Sekant:   G_s/G₀ = 1 / (1 + a|γ/γ₀.₇|),   a = 0.385   (Eq 7-3) → γ=γ₀.₇'de G_s=0.722 G₀
Kesme-gerilme: τ = G_s γ = G₀γ / (1 + 0.385 γ/γ₀.₇)            (Eq 7-7)
Tangent:  G_t = dτ/dγ = G₀ / (1 + 0.385 γ/γ₀.₇)²              (Eq 7-8)
```
γ = (3/2)ε_q (Eq 7-5); triaksiyalde γ = ε_axial − ε_lateral (Eq 7-6).

## 3. Alt-kesim (HS unload/reload rijitliğine)
```
G_t ≥ G_ur,  G_ur = E_ur/(2(1+ν_ur))                          (Eq 7-9)
γ_cutoff = (1/0.385)(√(G₀/G_ur) − 1) γ₀.₇                     (Eq 7-10)
```
γ ≥ γ_cutoff için G_t = G_ur (HS davranışına oturur). Yani HSsmall, küçük-şekildeğiştirmede E₀=2(1+ν_ur)G₀
ile başlar (≫ E_ur), γ_cutoff'ta E_ur'a iner → büyük şekildeğiştirmede HS ile AYNI.

## 4. Gerilme bağımlılığı + Masing
G₀ (ve dolayısıyla G_t), HS modülleriyle aynı güç yasasını izler: G₀ = G₀^ref·((c·cotφ−σ'₃)/(c·cotφ+p^ref))^m.
Masing (Eq 7-11): unload/reload'da γ₀.₇ → 2γ₀.₇ (eğri 2× rijit). HSsmall'da hardening plastisite virgin
loading'de hızlı degradasyonu zaten verir; unload/reload doğal olarak daha rijit (kısmen elastik).

## 5. FE entegrasyonu (Benz strain-history overlay)
HS'in elastik (unload/reload) rijitliği E_ur, şekildeğiştirmeye-bağlı E_t(γ_hist) ile DEĞİŞİR:
```
γ_hist = √3 ||H Δe|| / ||Δe||   (Eq 7-4; H = deviatorik şekildeğiştirme-tarihi tensörü, reversal'da resetlenir)
E_t(γ_hist) = 2(1+ν_ur) G_t(γ_hist)   (cut-off G_ur ile)
```
**Monotonik yükleme** (statik, reversal yok): H→I, γ_hist = birikmiş deviatorik kayma şekildeğiştirmesi.
GaussState'e γ_hist eklenir; her artımda Δγ kadar artar; elastik öngörücü E_t(γ_hist) kullanır. Tam cyclic
(H tensörü + Simpson-brick reversal) dinamik/cyclic için (dynamics ertelendi) — monotonik statik bu overlay'siz yeterli.

## 6. Doğrulama planı
1. **Degradasyon eğrisi (analitik):** G_s/G₀ ve G_t/G₀ vs γ/γ₀.₇ Eq 7-3/7-8'i birebir; γ₀.₇'de G_s=0.722G₀;
   γ_cutoff'ta G_t=G_ur; τ-γ hiperbolü Eq 7-7. (round-off)
2. **Tek-eleman:** HSsmall drained triaxial başlangıç rijitliği ≈ E₀ (≫ HS E₀=E_ur), γ_cutoff sonrası HS ile çakışır.
3. **Berlin kazı (MMM §17.7):** HS vs HSsmall yüzey-oturması farkı (HSsmall daha az/dar oturma çukuru).

İlgili: [[hardening-soil-formulation]], [[material-model-architecture]], [[plaxis-gap-analysis]] Faz A.5.
