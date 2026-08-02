# Mohr-Coulomb — Doğrulanmış Formülasyon (implementasyon spesifikasyonu)

Bu belge, KATAI 2D'nin MC bünye modelinde **uygulanacak tam denklemleri**, otoriter
kaynaklarla doğrulanmış biçimde sabitler. Amaç: kodlamadan önce matematiği kilitlemek
(felsefe: en basitten karmaşığa, her adım doğrulanır — yanlış gidersek sonradan
ayıklamak pahalı).

**Birincil kaynak:** Sysala & Čermák, *"Subdifferential-based implicit return-mapping
operators in Mohr-Coulomb plasticity"* (arXiv:1508.07435, 2016) — principal-stress
uzayında kapalı-form return; dört bölge + bölge karar kriterleri + tutarlı teğet.
**Destekleyici:** Clausen, Damkilde & Andersen (2007), *Computers & Structures* (aynı
kapalı-form yaklaşım); Adhikary et al. (2016), *IJNME* (genel çok-yüzeyli). PLAXIS-MM
(parametre/konvansiyon).

## 0. Konvansiyon (KRİTİK — solver'la tutarlı)
- **Tension-positive.** Basınç negatif (örn yüzey basıncı p → σyy=−p). Bütün solver
  bu konvansiyonda; MC de aynı olmalı.
- Asal gerilmeler **σ₁ ≥ σ₂ ≥ σ₃**; σ₁ en çekmeli, σ₃ en basınçlı.
- Mühendislik kayması γxy (B matrisi bu konvansiyonda).

## 1. Akma fonksiyonu ve plastik potansiyel
Sysala (3.2)–(3.3), bizim `mohr_coulomb.hpp` ile **birebir** (H=0, pekleşme yok):
```
f(σ) = (1+sinφ)·σ₁ − (1−sinφ)·σ₃ − 2c·cosφ            (akma; f≤0 kabul edilebilir)
g(σ) = (1+sinψ)·σ₁ − (1−sinψ)·σ₃                       (plastik potansiyel)
```
Apeks (tepe) hidrostatik çekme noktası: f=0 & σ₁=σ₂=σ₃ ⇒ **p_apex = c·cotφ**.

## 2. Elastik öngörücü (plane strain, εzz≡0)
λ, μ = Lamé; K = λ+2μ/3, G = μ. σzz açıkça taşınır (asal gerilmedir, τxz=τyz=0):
```
Δσxx = λ·(Δεxx+Δεyy) + 2μ·Δεxx
Δσyy = λ·(Δεxx+Δεyy) + 2μ·Δεyy
Δσxy = μ·Δγxy
Δσzz = λ·(Δεxx+Δεyy)
```
σ^tr admissible (f(σ^tr)≤0) ise adım elastiktir (return yok).

## 3. Return mapping — dört bölge (mükemmel plastik, H=0)
Notasyon: ν akış vektörü ∂g/∂σ; D_e:ν asal bileşenleri. Sysala (4.8)'den asal
güncelleme `σᵢ = σᵢ^tr − Δλ·[(2/3)(3K−2G)·sinψ + 2G·νᵢ]`. (2/3)(3K−2G)=2λ.

Bölge sınır değerleri (a priori, sadece σ^tr'den — Sysala 4.13/4.18/4.23):
```
γ_sl = (σ₁^tr − σ₂^tr) / (2G(1+sinψ))        γ_sr = (σ₂^tr − σ₃^tr) / (2G(1−sinψ))
γ_la = (σ₁^tr + σ₂^tr − 2σ₃^tr) / (2G(3−sinψ))
γ_ra = (2σ₁^tr − σ₂^tr − σ₃^tr) / (2G(3+sinψ))
```

### 3a. Düz yüzey (σ₁>σ₂>σ₃) — ν=(1+sinψ, 0, −(1−sinψ))
```
Δλ = f(σ^tr) / [ (4/3)(3K−2G)·sinψ·sinφ + 4G·(1+sinψ·sinφ) ]
σ₁ = σ₁^tr − Δλ[(2/3)(3K−2G)sinψ + 2G(1+sinψ)]
σ₂ = σ₂^tr − Δλ[(2/3)(3K−2G)sinψ]
σ₃ = σ₃^tr − Δλ[(2/3)(3K−2G)sinψ − 2G(1−sinψ)]
```
Geçerli ⇔ Δλ ≤ min(γ_sl, γ_sr).

### 3b. Sol kenar (σ₁=σ₂>σ₃; triaksiyel uzama) — ν₁+ν₂=1+sinψ, ν₃=−(1−sinψ)
```
num = ½(1+sinφ)(σ₁^tr+σ₂^tr) − (1−sinφ)σ₃^tr − 2c·cosφ
den = (4/3)(3K−2G)sinψ·sinφ + G(1+sinψ)(1+sinφ) + 2G(1−sinψ)(1−sinφ)
Δλ  = num/den
σ₁=σ₂ = ½(σ₁^tr+σ₂^tr) − Δλ[(2/3)(3K−2G)sinψ + G(1+sinψ)]
σ₃    = σ₃^tr − Δλ[(2/3)(3K−2G)sinψ − 2G(1−sinψ)]
```
Geçerli ⇔ γ_sl ≤ Δλ < γ_la.

### 3c. Sağ kenar (σ₁>σ₂=σ₃; triaksiyel basınç) — ν₁=1+sinψ, ν₂+ν₃=−(1−sinψ)
```
num = (1+sinφ)σ₁^tr − ½(1−sinφ)(σ₂^tr+σ₃^tr) − 2c·cosφ
den = (4/3)(3K−2G)sinψ·sinφ + 2G(1+sinψ)(1+sinφ) + G(1−sinψ)(1−sinφ)
Δλ  = num/den
σ₁    = σ₁^tr − Δλ[(2/3)(3K−2G)sinψ + 2G(1+sinψ)]
σ₂=σ₃ = ½(σ₂^tr+σ₃^tr) − Δλ[(2/3)(3K−2G)sinψ − G(1−sinψ)]
```
Geçerli ⇔ γ_sr ≤ Δλ < γ_ra.

### 3d. Apeks (σ₁=σ₂=σ₃=p)
```
p^tr = ⅓(σ₁^tr+σ₂^tr+σ₃^tr)
Δλ = [ (2/3)(σ₁^tr+σ₂^tr+σ₃^tr)·sinφ − 2c·cosφ ] / (4K·sinψ·sinφ)
p  = p^tr − 2K·sinψ·Δλ
```
Geçerli ⇔ Δλ ≥ max(γ_la, γ_ra). **Dejenere durum:** ψ=0 (veya φ=0) ⇒ payda 0;
akış deviatorik olduğundan tepe noktasına dönüş için doğrudan **p = c·cotφ** projeksiyonu
kullanılır (özel-durum).

### Bölge seçimi (Sysala Teorem 4.2 — sıralı kontrol)
1. `q_s(min(γ_sl,γ_sr)) < 0` ⇔ Δλ_smooth ≤ min(γ_sl,γ_sr) → **düz yüzey**.
2. değilse γ_sl<γ_la iken sol-kenar aralığı sağlanıyorsa → **sol kenar**.
3. değilse γ_sr<γ_ra iken sağ-kenar aralığı sağlanıyorsa → **sağ kenar**.
4. değilse → **apeks**.
(Pratikte: her bölgenin kapalı-form Δλ'sını hesapla, geçerlilik aralığını sağlayanı seç;
H=0'da q'lar afin olduğundan eşdeğer ve iterasyonsuz.)

## 4. Asal → tensör geri kurulumu (coaxiality)
σ ve σ^tr **aynı özvektörlere** sahiptir (Sysala 4.7 ardı). Plane strain'de σzz öz-yönü
ẑ; in-plane iki asal, in-plane öz-yönlerde. 3-yönlü sıralamada hangi asalın nereden
geldiği (in-plane-a / in-plane-b / zz) etiketlenir; return sonrası geri dağıtılır.
In-plane tensör, saklanan (cos2θ, sin2θ) yönelimiyle kurulur:
```
σxx = m + r·cos2θ,  σyy = m − r·cos2θ,  σxy = r·sin2θ      (m=(pa+pb)/2, r=(pa−pb)/2)
```

## 5. Tutarlı (algoritmik) teğet — P1.2c
Sysala Bölüm 5 (5.4–5.11): eigenprojeksiyon E_i ve türevleri ile bölge bazında 𝒟S;
çoklu özdeğer civarında dallanmasız (5.5/5.7/5.9/5.11). Plane strain sadeleştirmesi
Appendix A. Non-asosiye (ψ≠φ) ⇒ teğet **simetrik değil** → global Newton'da nonsimetrik
PARDISO kolu. (P1.2b'de geçici olarak elastik/continuum teğet; tutarlı teğet P1.2c.)

## 6. Doğrulama planı (Seviye 0 → 2)
- **Material-point (P1.2b):** dayatılan gerilme/şekildeğiştirme yolu; return sonrası
  f≈0 (yüzeye tam dönüş), analitik MC dayanımıyla (örn σ₁/σ₃ = Ka veya q_f) karşılaştırma,
  elastik round-trip, coaxiality.
- **P1.2c:** tutarlı teğet — sonlu farkla dσ/dε ≈ D_T sayısal doğrulaması.
- **Global (P1.3):** Prandtl Nc=5.14 (şerit temel) < %5.

## 7. Çekme kesmesi (tension cut-off) — Rankine kapağı (2026-07 eklentisi)

**Birincil kaynak:** PLAXIS MMM §3.2 Denk. (3-11) + §3.3.10: kesme açıkken üç ek akma
fonksiyonu `f_t,i = σᵢ − σ_t ≤ 0` (i=1,2,3) tanımlanır ve bu yüzeyler için **asosiye**
akış kullanılır; σ_t varsayılanı 0'dır ve MC modelinde kesme **varsayılan açıktır**.
Sıralı asallarda (σ₁ ≥ σ₂ ≥ σ₃) üç düzlem tek koşula iner: **σ₁ ≤ σ_t**.

**Kırpma:** MC kabul kümesinde en çekmeli asal apekste `c·cotφ`'dir (f≤0 & σ₁=σ₂=σ₃ ⇒
σ₁ ≤ c·cotφ) ⇒ `σ_t,eff = min(σ_t, c·cotφ)` (φ>0; φ=0'da kırpma yok). σ_t ≥ c·cotφ
girildiğinde kapak apekste teğettir ve hiç bağlanmaz (PLAXIS de girdiyi böyle sınırlar).
Safety/EC7 faktörlemesi c'yi böldüğünden kırpma otomatik daralır; σ_t'nin kendisi de
SRF/γ_c ile bölünür (MMM §4.3.7 Hoek-Brown'da Safety'nin cut-off değerini böldüğünü
açıkça yazar; MC'ye aynı kural — güvenli yön).

Birleşik kabul kümesi = MC piramidi ∩ {σ₁ ≤ σ_t}: dışbükey çokyüzlü. Tüm yüzeyler
asal uzayda DÜZLEM, akış yönleri sabit ⇒ her bölgede dönüş **afin**, J sabit-kesin,
teğet yine `principal_consistent_tangent` ile kurulur (FD gerekmez).

Notasyon: `D_e v = λ(Σv)·1 + 2G·v` (asal uzay); `D_e e₁ = (λ+2G, λ, λ)`;
`D_e m_mc = (2λsψ+2G(1+sψ), 2λsψ, 2λsψ−2G(1−sψ))` (m_mc = ∂g/∂σ, §3a ile aynı).
`s_m ≡ (σ_t(1+sφ) − 2c·cφ)/(1−sφ)` (MC yüzeyinin σ₁=σ_t düzlemiyle kesişimindeki σ₃;
kırpma gereği s_m ≤ σ_t).

### 7a. T yüzü (yalnız f_t1 aktif; asosiye n=e₁)
```
Δλ_t = (σ₁^tr − σ_t)/(λ+2G)
σ = (σ_t,  σ₂^tr − λΔλ_t,  σ₃^tr − λΔλ_t)
J: satır1 = 0;  satır2 = (−λ/(λ+2G), 1, 0);  satır3 = (−λ/(λ+2G), 0, 1)
```
Geçerli ⇔ σ₂^ret ≤ σ_t ve f_mc(σ^ret) ≤ 0. (f_mc değişimi = −Δλ_t[2λsφ+2G(1+sφ)] < 0
⇒ MC-kabul edilebilir trial'dan T-yüzü dönüşü MC'yi asla bozmaz.)

### 7b. TT kenarı (f_t1, f_t2; σ₁=σ₂=σ_t)
```
S = a+b = (σ₁^tr + σ₂^tr − 2σ_t)/(2λ+2G);  a−b = (σ₁^tr−σ₂^tr)/(2G)
σ = (σ_t, σ_t, σ₃^tr − λS)
J: satır1 = satır2 = 0;  satır3 = (−λ/(2λ+2G), −λ/(2λ+2G), 1)
```
Geçerli ⇔ b ≥ 0 (yüzey yetmedi), σ₃^ret ≤ σ_t, f_mc(σ_t, σ₃^ret) ≤ 0.

### 7c. T apeksi (f_t1..3; hidrostatik çekme köşesi)
```
σ = (σ_t, σ_t, σ_t);   J = 0
```
f_mc(σ_t,σ_t) = 2σ_t·sφ − 2c·cφ ≤ 0 kırpma gereği daima sağlanır.

### 7d. MC∩T doğrusu (düz MC yüzeyi + f_t1; non-asosiye MC akışı korunur)
Doğru = {(σ_t, s₂, s_m) : s_m ≤ s₂ ≤ σ_t} (V_R'den V_L'ye). Çarpanlar (a=Δλ_mc, b=Δλ_t):
```
A11 = 4λ·sψsφ + 4G(1+sψsφ)   (= §3a'daki payda)     A12 = 2λsφ + 2G(1+sφ)
A21 = 2λsψ + 2G(1+sψ)                                A22 = λ + 2G
[A11 A12; A21 A22]·[a; b] = [f_mc(σ^tr); σ₁^tr − σ_t]
σ₁ = σ_t,  σ₃ = s_m (kapalı-form),  σ₂ = σ₂^tr − a·2λsψ − b·λ
J: satır1 = satır3 = 0;  satır2 = e₂ᵀ − (2λsψ)·∂a/∂s − λ·∂b/∂s,
   [∂a/∂s; ∂b/∂s] = A⁻¹·[(1+sφ, 0, −(1−sφ)); (1, 0, 0)]
```
Geçerli ⇔ s_m ≤ σ₂^ret ≤ σ_t (ve a,b ≥ 0).

### 7e. Köşeler (nokta dönüşleri; J = 0)
```
V_L = (σ_t, σ_t, s_m)   [iki Rankine + iki MC yüzü — sol-kenar simetrisiyle 4 düzlem]
V_R = (σ_t, s_m, s_m)   [MC sağ kenarı ∩ T yüzü]
```

### 7f. Bölge seçimi — geçerlilik-kademeli (SS Koiter kenar deseni; toleranssız harita)
1. f_mc(σ^tr) ≤ tol ve σ₁^tr ≤ σ_t + tol → **elastik**.
2. f_mc(σ^tr) > tol → **MC-yalnız dönüş** (§3, DEĞİŞMEDEN — kesme kapalıyken bit-birebir
   eski yol). Kesme açık ve σ₁^ret ≤ σ_t + tol → kabul.
3. Aksi hâlde **T-yalnız kademesi** (7a → σ₂ aşarsa 7b → σ₃ aşarsa 7c); sonuç
   f_mc ≤ tol ise kabul. (Adım 2'den düşen her durumda σ₁^tr > σ_t'dir: MC dönüşü σ₁'i
   büyütmez ⇒ Δλ_t > 0 garantili.)
4. **Birleşik**: 7d doğrusu; σ₂^ret > σ_t → V_L; σ₂^ret < s_m → V_R.

Doğruluk: bölgeler trial uzayını döşer (mükemmel plastisite + lineer yüzeyler); kabul
koşulu = KKT (uygunluk + çarpan işaretleri) ⇒ ilk geçerli bölge doğru bölgedir; komşu
bölgeler ortak sınırda aynı noktayı verdiğinden harita süreklidir (Newton dostu).

### 7g. Doğrulama (test_mohr_coulomb + BVP)
- **Tersinme (her bölge):** bölge üzerinde nokta x* seç, trial = x* + Σλᵢ·D_e mᵢ (λᵢ>0,
  bölgenin aktif akış yönleri) kur ⇒ dönüş x*'ı birebir geri vermeli.
- **Kapalı formlar:** tek-eksenli uzama plane-strain'de akma sonrası σxx = σ_t,
  σyy = σzz = λσ_t/(λ+2G) DONUK; eş-iki-eksenli uzama → TT kenarı; hidrostatik → T apeks.
- **Teğet:** bölge içinde D_T vs merkezî FD; kesme-kapalı bit-birebir; BVP'de dönüş
  sonrası tarama σ₁ ≤ σ_t + tol.

Kaynaklar: arXiv:1508.07435 ; Clausen et al. 2007 (Comput. Struct.) ;
Adhikary et al. 2016 (IJNME, doi:10.1002/nme.5284) ; PLAXIS MMM §3.2 (3-11), §3.3.10.
