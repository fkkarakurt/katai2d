# Dinamik / sismik analiz — formülasyon (Faz D1: çekirdek + 1D serbest-alan)

KATAI 2D'nin dinamik/sismik izinin matematiksel temeli. Amaç: **birebir PLAXIS 2D Dynamics**
davranışı, kapalı-form 1D dalga teorisiyle < %5 doğrulanmış. Bu belge **çekirdeği** (hareket
denklemi + Newmark zaman integrasyonu + Rayleigh sönüm + sismik taban-uyarımı) ve **1D serbest-alan
site-response** doğrulama oracle'ını kilitler. 2D FEM kütle matrisi, absorbing/free-field sınırlar ve
TBDY entegrasyonu sonraki fazlar. İlgili: [[project-plaxis-parity-roadmap]].

**Kaynaklar (kilitli):** Newmark (1959) *J. Eng. Mech. Div. ASCE* 85:67; Chopra, *Dynamics of
Structures* (hareket denklemi, Newmark, Rayleigh, kararlı-hal); Kramer, *Geotechnical Earthquake
Engineering* (1996) §7 (1D dalga yayılımı, transfer fonksiyonu, site response); Hughes, *The Finite
Element Method* (tutarlı kütle, zaman integrasyonu); PLAXIS 2D 2025.1 Scientific Manual (Dynamics).
Konvansiyon: tension-pozitif, tüm kod tabanıyla aynı. TBDY 2018 sonraki faz.

## 1. Yarı-ayrık hareket denklemi (Chopra §9)
Uzayda FE ayrıklaştırılmış sistem:
```
M ü(t) + C u̇(t) + K u(t) = F(t)
```
- **M** = kütle matrisi. **Tutarlı (consistent):** Mₑ = ∫_Ω ρ NᵀN dΩ (eleman); ρ = kütle yoğunluğu
  = γ/g (γ birim ağırlık, g=9.81 m/s²). **Lumped (diyagonal):** satır-toplamı — açık dinamikte yaygın;
  implicit Newmark için tutarlı kütle daha doğru frekans verir (PLAXIS tutarlı kullanır). v1: tutarlı.
- **K** = statik rijitlik (mevcut `assemble_stiffness`; sismikte küçük-şekil-değiştirme elastik veya
  elastoplastik teğet). **C** = sönüm (bkz §3).
- **F(t)** = dış kuvvet; sismik taban-uyarımında efektif deprem kuvveti (bkz §4).

## 2. Newmark-β zaman integrasyonu (Newmark 1959; Chopra §15.2)
t_{n+1}=t_n+Δt için:
```
u_{n+1} = u_n + Δt·v_n + Δt²[(½−β)a_n + β·a_{n+1}]
v_{n+1} = v_n + Δt[(1−γ)a_n + γ·a_{n+1}]
```
γ=½, β=¼ (**ortalama-ivme / trapez kuralı**): koşulsuz kararlı, 2. mertebe doğru, **algoritmik sönüm
YOK** (enerji korur). (γ>½ sayısal sönüm ekler — yüksek-frekans gürültüsü bastırma; v1 γ=½.)
**Efektif rijitlik** (sabit Δt → BİR KEZ faktörle):
```
K_eff = K + a0·M + a1·C,   a0 = 1/(βΔt²),  a1 = γ/(βΔt)
```
**Efektif yük** (her adım):
```
F_eff = F_{n+1} + M(a0·u_n + a2·v_n + a3·a_n) + C(a1·u_n + a4·v_n + a5·a_n)
a2 = 1/(βΔt),  a3 = 1/(2β) − 1,  a4 = γ/β − 1,  a5 = Δt(γ/(2β) − 1)
```
K_eff u_{n+1} = F_eff çöz → sonra:
```
a_{n+1} = a0(u_{n+1}−u_n) − a2·v_n − a3·a_n
v_{n+1} = v_n + Δt[(1−γ)a_n + γ·a_{n+1}]
```
Başlangıç ivmesi: M·a0 = F(0) − C·v0 − K·u0. K_eff simetrik (K,M,C simetrik) → SPD ise PARDISO
mtype=2 / dense LU. **factor-once-solve-many** (konsolidasyon deseniyle aynı `SolveFactory` callback,
MKL-agnostik; boşsa dense LU referans yol).

## 3. Rayleigh sönüm (Chopra §11.4)
```
C = α·M + β·K
```
Modal sönüm oranı: ξ_n = α/(2ω_n) + β·ω_n/2  (ω_n = 2π f_n). İki hedef (f₁,ξ₁),(f₂,ξ₂)'den:
```
[1/(2ω1)  ω1/2] [α]   [ξ1]                       α = 2 ω1 ω2 (ξ1 ω2 − ξ2 ω1)/(ω2²−ω1²)
[1/(2ω2)  ω2/2] [β] = [ξ2]   →                    β = 2 (ξ2 ω2 − ξ1 ω1)/(ω2²−ω1²)
```
Eşit hedef ξ₁=ξ₂=ξ: α=ξ·2ω1ω2/(ω1+ω2), β=ξ·2/(ω1+ω2). Genelde iki hedef zemin tabakasının 1. ve 3.
doğal frekansı (ör. f₁ ve ~5f₁) seçilir → aradaki bandda ξ ≈ hedef, dışında artar (Rayleigh
sınırlaması; PLAXIS de aynı). `rayleigh_from_modes` bunu çözer.

## 4. Sismik taban-uyarımı (Kramer §7.2; Chopra §9.7) — relatif-deplasman formülasyonu
Rijit taban tek-eksenli a_g(t) ivmesiyle sarsılır. Toplam deplasman u_t = u + r·u_g (u = tabana GÖRE
relatif; r = etki vektörü: sarsma yönündeki serbest DOF'larda 1, diğerlerinde 0). u_t'yi denkleme
sokup rijit-cisim ötelemesinin K r u_g = 0, C r u̇_g ≈ 0 (öteleme rijit-mod → iç kuvvet üretmez)
olduğunu kullanınca:
```
M ü + C u̇ + K u = −M r a_g(t)      (efektif deprem kuvveti F_eff,eq = −M r a_g)
```
Yüzey **toplam** ivmesi ü_t = ü + a_g. Taban DOF'u sabitlenir (relatif u=0), F = −M r a_g tüm serbest
DOF'lara. (Esnek taban / uyumlu-taban dashpot'u ve absorbing sınırlar sonraki faz.)

## 5. 1D serbest-alan site-response — doğrulama oracle'ı (Kramer §7.2, altın standart)
Homojen elastik tabaka: kalınlık H, kayma dalga hızı Vs=√(G/ρ), taban rijit, yüzey serbest (SH dalga
düşey yayılır). Kayma kolonu = 1D çubuk (G↔E): ρ ü = G u,zz. FE: M ü + K u = f, K = ∫ G (dN/dz)²
(kayma), M = ∫ ρ N² (tutarlı; 2-düğüm lineer eleman → ρAL/6·[[2,1],[1,2]], k=GA/L·[[1,−1],[−1,1]]).
- **Doğal frekanslar (özdeğer K φ = ω² M φ):**  f_n = (2n−1)·Vs/(4H),  n=1,2,…  (temel f₁=Vs/4H,
  çeyrek-dalga rezonansı). M ve K montajını KESİN doğrular (dış çözücü gerekmez).
- **Sönümsüz transfer fonksiyonu** (yüzey/taban genlik oranı):  |F(ω)| = 1/|cos(ωH/Vs)|  →  f_n'de ∞.
- **Sönümlü** (frekans-bağımsız ξ):  |F(ω)| = 1/√(cos²(ωH/Vs) + (ξ·ωH/Vs)²).  Temel rezonansta
  (n=1) tepe **F_max ≈ 2/(π ξ)** (Kramer Eq 7.30). Bu, Newmark + Rayleigh + taban-uyarımı zincirini
  UÇTAN-UCA doğrular.

## 6. KATAI çekirdek uygulaması (Faz D1+D2 — test_dynamics)
- `analysis/dynamics.hpp`: `rayleigh_from_modes(f1,ξ1,f2,ξ2)` (§3); `solve_newmark(M,C,K,force(step),
  dt,nsteps,u0,v0[,γ,β,factory])` (§2, boyut-bağımsız → 2D için yeniden kullanılır). `NewmarkResult`:
  her adım u/v/a + t.
- **D1 doğrulama `test_dynamics` (a–c):** (a) **SDOF** (m,k,c): serbest-titreşim periyodu T=2π√(m/k) +
  sönümlü decay zarfı e^(−ξω t) + kararlı-hal harmonik dinamik-amplifikasyon faktörü Rd(β_f,ξ) —
  Newmark'ı KESİN kapalı-formla; (b) **1D kayma kolonu** özdeğer f_n=(2n−1)Vs/4H (M,K); (c) **taban-
  harmonik transfer fonksiyonu** |ü_t,yüzey/a_g| temel frekansta ≈ 2/(πξ) + f₁ konumu. Tümü ~KESİN
  (SDOF %0.00, kolon f₁ %0.01, transfer %0.0).

## 7. D2 — 2D FE tutarlı kütle (test_dynamics d/e)
- `assembly/assembler.hpp` `assemble_mass(mesh,dofs,density,builder)`: eleman-generic (tri6/tri15)
  tutarlı kütle Mₑ=∫ρNᵀN dA; skaler m_ij=ρ∫N_iN_j x VE y DOF bloklarına (blok-diyagonal, x-y kuplajsız).
  density=ρ=γ/g malzeme-başı. **Toplam kütle KESİN korunur** (1ᵀM1=2ρ·alan, ΣN_i=1 partition-of-unity;
  test (d): tri6/tri15 hata ~1e-14). Plane-strain; axisym kütle (r-ağırlıklı ∫ρNᵀN·r) sonraki faz.
- **2D serbest-alan doğrulama test (e):** SH kayma kolonu = uy≡0 (saf SH shear → σxx=0 kendiliğinden),
  rijit taban ux=0, yan+üst serbest, yatay taban-uyarımı F=−M·r·a_g. Kesin çözüm 1D kayma alanı ux=ux(y)
  → montajlı 2D M/K/C üzerinde AYNI `solve_newmark` 1D site-response'u üretir: özdeğer f_n=(2n−1)Vs/4H
  (%0.00) + transfer fn sub-rezonans 1/cos %0.1 + rezonans 2/(πξ) %0.0. **Boyut-bağımsız integratör 2D'de
  kanıtlandı** (aynı çekirdek 1D ve 2D FE'yi sürüyor).
- **Sonraki:** D3 (absorbing sınırlar, aşağıda §8), D4 (TBDY 2018 + GUI).

## 8. D3 — Absorbing (Lysmer-Kuhlemeyer viskoz) sınırlar (test_dynamics f/g)
Sonlu FE domaini bir yarım-uzayı temsil eder; dışarı giden dalga yapay sınırdan YANSIMAMALI (aksi
halde sahte enerji birikir). **Lysmer & Kuhlemeyer (1969) viskoz dashpot:** sınırda hıza-orantılı
çekme uygula:
```
t_damp = −[ c_n·v_n·n + c_t·v_t·t ],   c_n = ρ·Vp,  c_t = ρ·Vs
```
v_n=v·n, v_t=v·t (sınır normal/teğet hız); Vp=√((λ+2μ)/ρ) P-dalga, Vs=√(μ/ρ) S-dalga hızı. Sönüm
matrisine katkı: **C_b = ∫_Γ Nᵀ D_c N dΓ**, D_c = c_n(n⊗n) + c_t(t⊗t) (global 2×2, normal işaretinden
BAĞIMSIZ). `assemble_boundary_dashpot` (assembler) bunu kenar-integraliyle (4-nokta 1B Gauss, edge_shape)
toplar → global C'ye eklenir; sonra `solve_newmark` sönümü otomatik işler.
- **Yansıma katsayısı (1D normal geliş, KESİN):** empedans Z=ρ·V·A; sınır dashpot c için
  **R = (Z − c)/(Z + c)**. c=Z (Lysmer eşleşmesi) → R=0 (mükemmel yutma); c=0 (serbest) → R=+1; c=∞
  (sabit) → R=−1; c=Z/2 → R=1/3; c=2Z → R=−1/3. Yansıyan ENERJİ oranı = R².
- **D3 doğrulama `test_dynamics` (f):** 1D çubuk, d'Alembert TEK-YÖN aşağı pulse (u0=Gauss, v0=Vs·u0′),
  taban nodal dashpot c, tek-transit sonrası enerji E/E0 = R²: c=0→1.00, c=Z/2→1/9, c=Z→~0, c=2Z→1/9
  (yansıma formülü nicel doğrulandı). **(g):** 2D SH kayma kolonu (uy≡0), taban `assemble_boundary_dashpot`
  (c_t=ρVs → yatay SH hareketi teğet) → aşağı SH pulse tabandan RADYASYONLA çıkar (E_final/E0 küçük); sabit
  taban kontrolünde enerji hapsolur (korunur). Kenar-integral montajı kanıtlandı.
- **Sonraki:** D3b (free-field yan sınırlar, aşağıda §9), D4 (TBDY 2018 + GUI).

## 9. D3b — Free-field YAN sınırlar (test_dynamics h)
Yalnız absorbing dashpot yan sınırlar sismik girdide YANLIŞ olur: yandaki zemin sonsuza uzanır ve
serbest-alan (1D site-response) hareketini yapar — sınır bunu SÜRDÜRMELİ, sıfıra yutulmamalı. Aynı
zamanda içeriden saçılan (yapı/eğik-geliş) dalgaları yutmalı. **Lysmer free-field (PLAXIS/FLAC):**
yan sınıra serbest-alan kolonu (u_ff) dashpot ile bağlanır; ana-ızgaraya uygulanan net kuvvet
```
F = C_b·(v_ff − v_2D)          (SH: yalnız teğet; σxx^ff = 0 çünkü u_ff x'ten bağımsız → iç kuvvet yok)
```
İki parçaya AYRIŞIR: **−C_b·v_2D → dashpot** (global C'ye, absorbing gibi) + **+C_b·v_ff → sürücü
kuvvet** (her adım force(step)'e; v_ff serbest-alan sınır hızı). Yeni çekirdek gerekmez: D3'ün
`assemble_boundary_dashpot`'u (C_b) + solve_newmark force callback'i. Serbest-alan hareketindeyken
(v_2D=v_ff) net kuvvet SIFIR → site-response bozulmaz; saçılma (v_2D≠v_ff) → dashpot yutar.
- **D3b doğrulama `test_dynamics` (h):** 2D SH domain (uy≡0), rijit taban + Ricker pulse a_g(t). **A**=serbest
  yan (referans serbest-alan, yanal-üniform → u_A=u_ff). **B**=free-field yan (C=C_b, force += C_b·v_A) →
  u_B ≈ u_A (KESİN; serbest-alan korunur, sürücü+dashpot iptal olur). **C**=yalnız absorbing yan (C_b, sürücü
  YOK) → u_C serbest-alan hareketini YANLIŞ yutar, yanlarda u_A'dan belirgin sapar. B≪C → free-field sürücü
  kuvvetinin şart olduğu kanıtlandı. (Lineerlik: sürücü serbest-alanı korur + dashpot sapmayı §8'de yutar →
  genel durumu kapsar.)
- **Sonraki:** D4 = TBDY 2018 (Türk tasarım spektrumu S_DS/S_D1, zemin sınıfları ZA–ZE, spektral eşleştirme,
  sismik tasarım kontrolleri) + GUI (Dynamic faz tipi, ivme-kaydı girişi, response-spectrum). NOT: D4 ivme
  zaman-serisi + ara-adım saklama → ertelenen P3 zaman-playback ortak altyapı (D2 her adım u/v/a saklıyor →
  2D'de selektif kayıt gerekebilir).

## 10. D5 — Zemin-yapı etkileşimi (SSI): yapı dinamik sistemde (`analysis/structural_dynamics.hpp`)
D4b'ye kadar dinamik faz ZEMİN-ONLY idi: yapısal elemanlar DofMap'e DOF ekliyor (plate dönme φ, bağımsız
perde ötelemesi) ve mesh'i bölüyordu, ama K/M'ye GİRMİYORDU → o DOF'ların satırları tamamen boş → K_eff
TEKİL → çözücü heap'i bozuyordu. D5 yapıyı sisteme kurar: zemin ve duvar BİRLİKTE sallanır.

### 10.1 Rijitlik — statik teğetin elastik dalı
`assemble_structural_stiffness` zemin K'sının ÜSTÜNE yapısal elastik rijitliği ekler (öteleme DOF'ları
PAYLAŞILDIĞI için montaj toplaması kuplajı kendiliğinden kurar):
```
plate (3/5 düğüm):  K_p = ∫(EA BεᵀBε + EI BκᵀBκ)ds + ∫kGA' BγᵀBγ ds   (seçmeli azaltılmış; §plate)
ankraj:             K_a = (EA/L)·(g⊗g),           g = ∂U/∂u
geogrid:            K_g = ∫ EA BᵀB ds             (2-nokta Gauss)
interface:          K_i = Σ_q w_q J_q [k_s(a⊗a) + k_n(b⊗b)]   (Newton-Cotes düğüm-çiftleri, Day&Potts)
```
Bu, `solve_nonlinear`'ın u=0'daki teğetinin ta kendisidir → **kimlik test edilir** (test_ssi_dynamics (a),
9.5e-15). **v1 LİNEER:** geogrid tension-only kesimi, ankraj akması, interface Coulomb kayması dinamik dalda
YOK — gergin/kaymayan durum etrafında lineerleştirme (sınırlar: `docs/validation/seismic-verification.md`).

### 10.2 Kütle — plate ataleti
Yalnız plate kütle taşır (ankraj/geogrid/interface PLAXIS'te de ağırlıksız):
```
M_öteleme = ∫ ρA NᵀN ds   (u_x ve u_y blokları ayrı, kuplajsız)      ρA = w/g      [Mg/m]
M_dönme   = ∫ ρI NᵀN ds   (φ DOF'u)                                  ρI = ρA d²/12 [Mg·m]
```
d = eşdeğer kalınlık √(12EI/EA) (PLAXIS Ref. §5.6: aynı EA ve EI'ye sahip dikdörtgen kesit) → ρI o kesitin
kütlesel ikinci momenti. Doğrulama: konsol f_n=(β_nL)²/(2πL²)√(EI/ρA) (Blevins Tablo 8-1) %0.01.

### 10.3 Etki vektörü r — rijit taban ötelemesi
```
r = 1  {zemin düğümü u_x} ∪ {BAĞIMSIZ perde öteleme u_x};   r = 0  {düşey}, {dönme φ}
```
Perde bağımsız DOF'u ATLANIRSA duvar tabandan atalet kuvveti ALMAZ → sessizce yanlış SSI. Sınanır:
**K·r = 0** (rijit öteleme sıfır-enerji modu — teorem, sayı değil) + **rᵀMr = toplam öteleme kütlesi**.
*Fixed-end ankraj:* uzak uç relatif çerçevede sabittir = tabanla birlikte hareket eder (zemine/kayaya
ankastre kök idealizasyonu). Mutlak çerçevede u_g SADELEŞİR → elongasyon yalnız relatif alandan gelir;
formülasyon tutarlıdır.

### 10.4 TEKİL KÜTLE — başlangıç ivmesi (SSI'ye özgü çökme sınıfı)
Yapı girince **M genelde TEKİLDİR**: dönme DOF'u yalnız ρI taşır, ve **w=0 (GUI varsayılanı) plate hiç
kütle katmaz** → M'de tam-sıfır satırlar. Newmark bunu KALDIRIR (M yalnız matvec; K_eff = K+a₀M+a₁C tekil
değil çünkü K her yapısal DOF'a rijitlik verir). Tek sorun VARSAYILAN başlangıç ivmesidir: M a = F(0)−Cv₀−Ku₀
→ M'yi FAKTÖRLER. Tekil M'yi PARDISO SPD'ye vermek erişim ihlali üretir (C++ exception DEĞİL → `/EHsc
catch(...)` yakalamaz → mesajsız kapanma). **Çözüm kapalı-form, yaklaşıklık değil:** dinlenmeden başlarken
```
a(0) = −r·a_g(0)   ⇒   M a(0) = −M·r·a_g(0) = F(0)   HER M için KESİN (tekil olsa bile)
```
Fiziği: t=0'da sistem henüz deforme olmamıştır → relatif çerçevede tümü rijit olarak −a_g(0) ile ivmelenir
(toplam ivme a_rel+a_g = 0). Free-field'de de geçerli: v_ff(0)=0 (kolon da dinlenmeden başlar). Verilince M
hiç faktörlenmez (`solve_newmark` `a0_init` parametresi; tekil-olmayan M'de faktörleyen yolla 9.9e-17 örtüşür).


## 11. COMPLIANT BASE (absorbing taban) — KİLİTLİ FORMÜLASYON (v0.6; İMPLEMENTASYON BEKLİYOR)

**Birincil kaynak, birebir:** PLAXIS 2D 2025.1 Scientific Manual §6.3.1-6.3.2 + Tutorial §17.8.5.

### 11.1 Kaynaktan kilitlenen konvansiyonlar
- **Viskoz sınır (Lysmer & Kuhlemeyer 1969, Sci §6.3.1):** sınırda sabitlik yerine damper;
  normal/kayma bileşenleri c_p = ρV_p, c_s = ρV_s (birim alan başına), gevşeme katsayıları C1=C2=1
  "pratik uygulamalar için yeterli".
- **Compliant base (Joyner & Chen 1975, Sci §6.3.2, Eqn 6-12):** yutma ve dinamik girdi AYNI yerde,
  taban sınırında. Gerilme, YUKARI ve AŞAĞI parçacık hızlarıyla yazılır ve — birebir alıntı —
  *"The reaction of the dashpots is multiplied by a factor 2 since half of the input is absorbed by
  the viscous dashpots and half is transferred to the main domain."* C1=C2=1 ŞART ("works correctly
  if the relaxation coefficients are equal to 1"). Yani taban kayma traksiyonu:
  τ(t) = ρ V_s [ 2 u̇↑(t) − u̇(t) ]   (u̇↑ = girdi yukarı-giden dalga; u̇ = ana domain toplam hızı).
- **Girdi yarılama (Tut §17.8.5, birebir):** *"Considering that the boundary condition at the base of
  the model will be defined using a Compliant base, the input signal has to be taken as half of the
  bedrock (within) motion"* (ux,start,ref = 0.5 ile). Yorum: kayıt taban seviyesindeki TOPLAM (within)
  hareket kabul edilir; yukarı-giden bileşen = yarısı. Rijit-taban yolundaki "within = girdi aynen"
  konvansiyonundan FARKI ürün içinde açıkça yazılacak (sessiz yarım/iki-kat hatası = klasik tuzak).

### 11.2 KATAI implementasyon planı (kilitli)
1. Opt-in `Phase.seismic_compliant_base` (varsayılan KAPALI = rijit taban, bit-birebir korunur).
2. Formülasyon TOPLAM harekette: taban x-DOF'ları serbest; C'ye taban dashpotları (ρV_s, tributary
   uzunlukla) eklenir; F(t) taban düğümlerinde 2ρV_s·u̇↑(t) traksiyonu. Rijit-taban yolundaki
   −M·r·a_g relatif sürücüsü bu dalda KULLANILMAZ (taban hareketi artık bilinmeyen).
3. Yanal free-field kolonları da compliant tabanlı çözülmeli (1D kolon + aynı taban koşulu), yoksa
   yan sınır rijit-taban serbest-alanını dayatır = tutarsızlık.
4. Monitör/spektrum: toplam formülasyonda a_yüzey doğrudan toplamdır (a_g eklenmez — rijit yol
   a_rel + a_g ekliyordu; karışırsa sessiz çift sayım).
5. V&V oracle'ı HAZIR DESEN: test_site_response_benchmark'ın katmanlı transfer-matris oracle'ı,
   taban BC'si u=0 yerine RADYASYON koşuluyla (G*û′ = iωρ_rV_s,r(2û↑ − û), tam Joyner-Chen biçimi)
   toplam-hareket olarak genişletilir. Sönümsüz tek tabakada kapalı-form pin:
   |T_outcrop(ω)| = 1/√(cos²(kH) + α² sin²(kH)),  α = ρ₁V_s,1/(ρ_r V_s,r)  → rezonansta |T| = 1/α
   (radyasyon sönümü; rijit tabanın sonsuz tepesinin yerine SONLU tepe — fiziksel fark budur ve
   testin dişi de budur: rijit-taban koşusu aynı profilde 1/α'yı AŞMALI, compliant koşu OTURMALI).