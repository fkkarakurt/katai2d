# Efektif Gerilme & Boşluk Basıncı — Formülasyon (uygulama spesifikasyonu)

Geoteknikteki en temel fizik: zemin davranışı **efektif gerilme** ile yönetilir.
Kaynaklar: Terzaghi (1923) efektif gerilme ilkesi; Biot (1941) konsolidasyon; PLAXIS
Reference/Material Models Manual (drained/undrained/consolidation). Felsefe: önce
matematiği kilitle (bkz `mohr-coulomb-formulation.md`), sonra en küçük parçayı doğrula.

## 0. Konvansiyon (solver'la tutarlı)
- **Tension-positive** (tüm solver böyle). Boşluk basıncı **u ≥ 0 bir basınçtır**
  (çekme-pozitif konvansiyonda gerilmeye −u katkı yapar).
- Terzaghi: efektif = total − boşluk basıncı katkısı. Çekme-pozitifte:
  ```
  σ_total = σ' − u·m,      yani   σ' = σ_total + u·m
  ```
  m = boşluk-basıncı yön vektörü: **plane strain** m = [1, 1, 0] (σxx, σyy normal;
  σxy kayma etkilenmez), **axisymmetric** m = [1, 1, 0, 1] (σr, σz, σθ normal).
- Sıkışmaz su varsayımı (Terzaghi): u tüm normal bileşenlere eşit etki eder.

## 1. Denge — efektif gerilme + boşluk basıncı yükü
Denge **total gerilmede** kurulur:
```
∫ Bᵀ σ_total dV = f_ext
∫ Bᵀ (σ' − u·m) dV = f_ext
∫ Bᵀ σ' dV = f_ext + ∫ Bᵀ u·m dV
```
**Anahtar gözlem (önceden-tanımlı u için):** bünye modeli σ' (efektif) üretir; boşluk
basıncı sağ-tarafa **ek bir yük** olarak girer:
```
f_pore = ∫ Bᵀ u·m dV       (her Gauss noktasında u·m, B ile)
```
Yani **malzeme ve solver DEĞİŞMEZ** — yalnızca yük vektörüne f_pore eklenir ve
hesaplanan "gerilme" artık **efektif gerilmedir** (çünkü yük boşluk terimini içerir).
Bünye (MC akma) doğal olarak σ'·ile çalışır. Gravity **doygun birim ağırlık γ_sat**
ile (total denge).

Plane strain'de düğüm i katkısı (m=[1,1,0], Bᵀm = [∂Nᵢ/∂x, ∂Nᵢ/∂y]):
```
f_pore[2i]   += w·detJ · u · ∂Nᵢ/∂x
f_pore[2i+1] += w·detJ · u · ∂Nᵢ/∂y
```
(Axisymmetric: m=[1,1,0,1], r-ağırlıklı, Bᵀm = [∂Nᵢ/∂r + Nᵢ/r, ∂Nᵢ/∂z].)

## 2. Boşluk basıncı alanı (ilk hedef: hidrostatik / steady-state)
Su tablası kotu z_w; hidrostatik:
```
u(x, z) = max(0, γ_w · (z_w − z))      (su tablası altında; üstünde 0, suction yok)
```
γ_w = suyun birim ağırlığı (~9.81 kN/m³). Bu, prescribed (deplasmandan bağımsız)
alandır → f_pore sabittir, teğeti etkilemez.

## 3. Doğrulama (Seviye 0 → 1)
- **Batık kolon (elastik):** su tablası üstte (z_w = H), doygun zemin, γ_sat gravity +
  hidrostatik u. Beklenen: **efektif düşey gerilme σ'_v = −γ'·(H−z)** (kaldırma kuvveti),
  γ' = γ_sat − γ_w (gömülü/buoyant birim ağırlık); boşluk basıncı u = γ_w·(H−z).
  Round-off (lineer alan, FE uzayında).
- (Sonra) **MC + su:** undrained değil, drained efektif-gerilme analizi; örn şev FoS
  su tablasıyla (efektif dayanım c', φ'·ile).

## 4. Undrained (A) — excess pore pressure üretimi (PLAXIS Undrained A)

Ani yükleme: su boşalamaz, near-incompressible davranır ve **excess** boşluk basıncı
üretir. PLAXIS "Undrained (A)": **efektif** stiffness (E', ν') + **efektif** dayanım
(c', φ') girilir; suyun bulk modülü iskelete eklenir.

**Su bulk modülü / porozite (varsayılan νu yaklaşımı, PLAXIS):**
```
K_w/n = ( 3(ν_u − ν') / ((1−2ν_u)(1+ν')) ) · K',     K' = E' / (3(1−2ν'))
```
ν_u = undrained Poisson (PLAXIS varsayılan **0.495**; tam 0.5 stiffness'i tekil yapar).
Türev (doğrulandı): bu, doğru undrained bulk modülünü verir:
```
K_u = K' + K_w/n = E'(1+ν_u) / (3(1+ν')(1−2ν_u))       (= 2G(1+ν_u)/(3(1−2ν_u)))
```

**Undrained (total) stiffness:** efektif operatöre su volumetrik rijitliği eklenir:
```
D_u = D'(E', ν') + (K_w/n) · m·mᵀ
```
m = [1,1,0] (plane strain) / [1,1,0,1] (axisym). m·mᵀ yalnız normal-normal bloku
doldurur (kayma etkilenmez).

**Excess boşluk basıncı (volumetrik şekildeğiştirme kuplajı):**
```
Δu_excess = −(K_w/n)·Δε_v ,   ε_v = mᵀε = ε_xx+ε_yy(+ε_zz)   (tension-pozitif: basınç ε_v<0 → u>0)
```
Toplam boşluk basıncı u = u_steady (§2) + u_excess. Total gerilme σ = σ' − u·m.
Bünye (MC) **efektif** σ'·ile çalışır; D_u yalnız global rijitlik/teğette kullanılır.

**Skempton B (doygunluk göstergesi):**
```
B = (K_w/n) / (K' + K_w/n) = (K_w/n) / K_u
```
ν_u=0.495, ν'=0.3 → K_w/n = 45·K' → **B = 45/46 ≈ 0.978** (PLAXIS varsayılanı; tam 1
için ν_u→0.5 ama tekillik). Doygun zeminde B≈1 beklenir.

**Uygulama notu (TAMAM):** Global undrained `solve_nonlinear`'a bağlandı. Bünye **efektif**
σ' ve efektif teğet üretir; solver malzeme `undrained` ise her Gauss noktasında **total**
gerilme σ = σ' + (K_w/n)ε_v·m kullanır, global teğete (K_w/n)m·mᵀ = D_u ekler ve excess
boşluk basıncını `GaussState::eps_vol` içinde taşır (u = −(K_w/n)ε_v). Aynı kod yolu hem
LinearElastic hem MC için çalışır (MC'de efektif σ' return + teğet augmentasyonu).
Doğrulama `test_undrained_global`: konfine (oedometer) undrained yükleme → 1D analitik
M_u=M'+K_w/n; **u/q=(K_w/n)/M_u=0.9653** (ν_u=0.495), σ'_yy=−M'q/M_u, total σ_yy=−q
denge, oturma=−qH/M_u (drained'e göre ~29× rijit) — round-off. (Klasik Terzaghi t=0
koşulu: ani yüklemede uygulanan gerilmenin neredeyse tamamını boşluk suyu taşır.)
- **Konsolidasyon (Biot):** eşlenik akış-deformasyon, zaman-bağımlı; boşluk basıncı
  ek DOF (u alanı). Büyük ek.
- **Staged construction'da su:** her fazda farklı su tablası / boşluk basıncı.

Kaynaklar: Terzaghi (1923); Biot (1941) *J. Appl. Phys.* 12:155–164; PLAXIS Reference &
Material Models Manual (Bentley); [[literature-review]].
