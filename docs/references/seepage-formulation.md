# Yeraltısuyu Akışı (Steady-State Seepage) — Formülasyon (P3 / desteklenmiş kazı için)

Kazı/perde problemlerinde su akışı (drawdown, kaldırma, sızma) boşluk basıncı dağılımını belirler;
bu dağılım efektif gerilme analizine girer (σ'=σ−u). PLAXIS'in "groundwater flow" adımı. Önce
**confined (sınırlı) steady-state** akış (lineer Laplace); sonra **unconfined** (serbest yüzey, nonlineer).

**Kaynaklar (kilitli):**
- **Bear (1972)**, *Dynamics of Fluids in Porous Media* — Darcy yasası, süreklilik.
- **Harr (1962)**, *Groundwater and Seepage* — klasik akış ağları (flow net), analitik çözümler.
- **Smith, Griffiths & Margetts**, *Programming the FEM* (5. baskı) — çalışan FE sızma kodu (Laplace/
  Poisson, confined + serbest yüzey). Doğrudan implementasyon referansı (SGM).
- **PLAXIS 2D Reference Manual** — groundwater flow (steady-state + transient), boundary tipleri.

## 0. Konvansiyon
2D düzlem akış (birim kalınlık). **Hidrolik yük** (head) h = y + u/γ_w (kot + basınç yükü); y = kot
(yukarı +), u = boşluk (basınç) basıncı, γ_w = suyun birim ağırlığı (~9.81 kN/m³). Boşluk basıncı:
```
u(x,y) = γ_w · (h − y)        (h > y, su tablası altında);   u = 0   (h ≤ y, suction yok)
```
Tension-pozitif solver'la tutarlı: u ≥ 0 basınç; efektif gerilmeye σ'=σ+u·m (bkz effective-stress).

## 1. Darcy + süreklilik (yöneten denklem)
Darcy özgül debisi (specific discharge): q = −k·∇h, k = geçirgenlik (permeability) tensörü
[[k_x,0],[0,k_y]] (ana eksende anizotropik; birim: uzunluk/zaman). Sıkışmaz, steady-state süreklilik:
```
∇·q = 0   ⇒   ∇·(k ∇h) = 0          (anizotropik Laplace; k sabitse k_x h_xx + k_y h_yy = 0)
```

## 2. Galerkin FE
Yük alanı h ≈ Σ N_i h_i (aynı tri6/tri15 şekil fonksiyonları). Zayıf form:
```
∫_Ω (∇w)ᵀ k (∇h) dΩ = ∮_Γ w (k∇h·n) dΓ = −∮_Γ w q_n dΓ
```
**Eleman geçirgenlik (conductivity) matrisi:**
```
H_e = ∫_Ω Gᵀ k G dΩ ,   G = [∂N/∂x ; ∂N/∂y]  (2×n; ∂N_i/∂x = B(0,2i), ∂N_i/∂y = B(1,2i+1))
```
G satırları mevcut `strain_displacement` B'sinden çıkarılır (eleman-generic). Global sistem:
```
H · h = Q
```
Q = düğümsel debi (inflow +). H simetrik pozitif tanımlı (k>0) → PARDISO SPD.

## 3. Sınır koşulları
- **Dirichlet (verilen yük):** h = h₀ (rezervuar/su seviyesi, sızma yüzeyi çıkışı). İnhomojen ⇒
  partition: serbest sistem H_ff h_f = Q_f − H_fp h_p ("Dirichlet lift"; sabit düğüm kolonu sağa).
- **Neumann (debi):** q_n verilen → Q_i += ∫ N_i q_n ds (içeri-akış pozitif), `assemble_seepage_flux`
  (kenar integrali, yüzey çekmesiyle aynı 1B Gauss). Kuyu/infiltrasyon/verilen-debi. **Geçirimsiz
  (impermeable) sınır q_n=0 = DOĞAL** (hiçbir şey yapma). *(test_seepage: 1D kolon h₀+sağ akı q_n →
  h(x)=h₀+(q_n/k)x round-off, sol çıkış = −q_n·H, korunum)*
- **Serbest yüzey (unconfined):** bilinmeyen freatik yüzey → nonlineer (bkz §6).
- **Seepage face (sızma yüzeyi):** serbest-drenajlı mansap yüzünde, çıkış noktası üstünde h=y (basınç=0);
  aktif-set iterasyonu (sonraki — desteklenmiş kazı drenajı için).

## 4. Deformasyonla kuplaj (efektif gerilme)
Akış çözümünün düğümsel head alanı → boşluk basıncı u(x,y)=γ_w·max(0,h−y) → boşluk yükü f_pore=∫Bᵀu·m dV
→ efektif gerilme analizi (malzeme/solver DEĞİŞMEZ, bkz effective-stress-formulation.md §1). Akış ve
deformasyon **ayrık (uncoupled)** çözülür (steady-state; konsolidasyon=Biot eşlenik, sonraki).
Üretim devri `assemble_pore_load_from_head` (seepage.hpp): düğümsel head'i her Gauss noktasında interpole
eder (ψ=h_gp−y_gp, u=γ_w·max(0,ψ)), konfine (ψ≥0) tam interpole. PLAXIS'in "groundwater flow → effective
stress" devrinin birebir karşılığı.

**Terzaghi yukarı-sızma (kuplaj doğrulaması, kapalı-form):** doygun kolon, sabit yukarı gradyan i,
lineer head h(y)=H+i(H−y). Efektif gerilme (Terzaghi 1943; Das, *Principles of Geotechnical Engineering*;
Craig, *Soil Mechanics*):
```
σ'_v(y) = −(γ' − i·γ_w)(H − y),   γ' = γ_sat − γ_w        (tension-poz; z=H−y derinlik)
```
Sızma kuvveti i·γ_w birim-derinlik başına (batık) efektif gerilmeyi azaltır. **Kritik gradyan
i_cr=γ'/γ_w'de σ'_v=0** = kaynama/heave (quick condition) — kazı taban stabilitesinin temel kriteri.
i=0 → −γ'(H−y) (batık, §1 ile tutarlı). Statik-belirli ⇒ E,ν'den bağımsız, FE round-off üretir.
*(test_seepage: i=0 / i=0.3 / i=i_cr üç rejim, max|σ'_v err|~1e-12)*

## 5. Doğrulama planı (basitten karmaşığa)
1. **1D Darcy kolonu (confined):** yatay/dikey kolon, iki uçta verilen yük h₁,h₂, yanlar geçirimsiz →
   **lineer yük profili** h(x)=h₁+(h₂−h₁)x/L (round-off), üniform debi q=k(h₁−h₂)/L. *(test_seepage)*
2. **2D confined patch:** lineer yük alanı h=a·x+b·y sınırda verilince içeride birebir üretilir
   (yama testi); anizotropik k. *(test_seepage)*
2b. **MMS yakınsama (eğri alan kesinliği):** patch test *yapısı gereği* lineer alanı birebir verir
   → çözücünün eğri bir alanda DOĞRU + teorik hızda yakınsadığını kanıtlamaz. **Manufactured
   Solution** (Roache 1998): **h = eˣ·sin(y)** HARMONİKtir (∇²h = eˣsiny − eˣsiny = 0), izotropik k
   için kesin sızma çözümü. **Transandantal** (polinom değil) seçildi: NE tri6 NE tri15 birebir
   üretebilir → ikisi de inceldikçe kendi teorik hızında yakınsar: **tri6 (P2) ~O(h³)** (gözlem
   3.56→3.79, nodal süperyakınsama), **tri15 (P4) ~O(h⁵)** (gözlem 4.68→4.87). tri15 8×4 (4.7e-7),
   tri6 24×12'den (1.4e-6) çok daha az elemanla daha doğru. **NOT:** harmonik *polinom* (örn
   x³−3xy² = Re((x+iy)³)) Laplace Galerkin çözücüsünde düğüm değerlerinde round-off'a yakalanır
   (kuadratik tri6 bile), bu yüzden yaklaşım hatasını SINAMAZ — transandantal alan şart. *(test_seepage)*
2c. **Radyal confined akış (ilk gerçek 2B eğri-alan + debi):** önceki testler 1B (Darcy/Terzaghi) ya
   da manufactured. **Thiem/well-flow** kapalı-form (Harr 1962; Bear 1972; Thiem 1906):
   ```
   h(r) = h₁ + (h₂−h₁)·ln(r/r₁)/ln(r₂/r₁) ,   Q = k·(h₁−h₂)·θ/ln(r₂/r₁)   (r'den bağımsız ⇒ korunum)
   ```
   = elastikteki Lamé silindirinin (test_axisym_cylinder) sızma karşılığı. **Mesh:** mantıksal
   (r,θ) dikdörtgeninde yapısal tri6 → düğümler fiziksel (r·cosθ, r·sinθ)'ye remap (izoparametrik
   eğri kenarlar; yeni mesher YOK). **Debi geri-kazanımı `compute_nodal_flux` (Q=K·h):** Dirichlet
   düğümlerinde sınır debisi, serbest düğümlerde ≈0; TÜM düğüm toplamı=0 (kütle korunumu, K satır
   toplamı sıfır). Sonuç: head order **2.91 (~O(h³))**; **debi relerr 1.28e-6** (consistent flux
   süperyakınsar — sınır debisi varyasyonel olarak çok daha doğru); kütle korunumu ~1e-18. *(test_seepage)*
3. **Confined akış ağı — baraj altı sızma + uplift (kesin arccos yasası):** düz tabanlı baraj,
   homojen izotrop yarı-sonsuz temel; taban boyunca yük **arccos yasası** (Harr 1962; Polubarinova-
   Kochina; Das):
   ```
   h(x) = (H/π)·arccos(x/b) ,   −b ≤ x ≤ b      (memba H @ −b, mansap 0 @ +b, merkez H/2)
   ```
   **Perde = üst sınırda geçirimsiz segment** (içsel slit GEREKMEZ — klasik dam-base kurulumu).
   **Kenar tekilliği:** exit gradient baraj topuğunda (x=±b) SONSUZ (fiziksel; her FE kodunda, PLAXIS
   dahil, mesh-duyarlı → mühendislikte integral/ortalama kullanılır). Doğrulama: pürüzsüz **merkezi
   taban** (|x|≤0.5b) arccos'a yakınsar (kaba %0.70 → ince %0.49); h(0)=H/2 ve uplift kuvveti ∫h dx=H·b
   ANTİSİMETRİden (h(x)+h(−x)=H) kesin (mesh-bağımsız) ⇒ arccos ŞEKLİ asıl test. Kütle korunumu ~1e-16.
   **Teşhis bulgusu:** sapma TRUNCATION değil (domain 12→24 büyütünce 0.72→0.69%) DISCRETIZATION/tekillik
   baskın (mesh 0.25→0.125 ⇒ 0.72→0.47%). *(test_seepage)*
   KALAN: sonlu-derinlik debi (Harr eliptik-integral form-faktörü) + sheet-pile penetrasyonu.
4. **Unconfined (freatik yüzey):** nonlineer iterasyon (freatik üstü k azaltma); kazı drawdown.

**Kaynak (V&V):** Roache, P.J. (1998), *Verification and Validation in Computational Science and
Engineering* — Method of Manufactured Solutions; gözlemlenen yakınsama mertebesinin teorik
değere (P2→3, P4→5) eşitliği çözücü doğruluğunun altın-standart kanıtıdır.

## 6. Unconfined (serbest-yüzeyli) akış — değişken bağıl geçirgenlik
Freatik yüzeyin (su tablası) domain içinde olduğu durum. **Sabit mesh, değişken k_rel(ψ)** yaklaşımı
(Bathe & Khoshgoftaar 1979; PLAXIS akış modülü relative-permeability): basınç yükü ψ=h−y; ψ≥0 doygun
(k_rel=1), ψ<0 doymamış (k_rel→k_min). Eleman iletkenliği He=∫Gᵀ k·k_rel(ψ_gp) G; freatik yüzey ψ=0
konturunda kendiliğinden oluşur. **Picard iterasyonu** (`solve_unconfined_seepage`): her adımda
k_rel(head) güncelle → SPD çöz → alt-gevşet (under-relax).

**k_rel modeli — van Genuchten (transient/coupled ile tutarlı):** `UnconfinedOptions.retention`
verilirse (GUI yolu: `build_flow.hpp` her malzemenin `gw_ga/gw_gn/gw_gl/gw_Sres`'inden kurar) k_rel =
**van Genuchten/Mualem** (`materials/water_retention.hpp`, `relative_permeability_psi`), yani transient
+ fully-coupled akışın kullandığı AYNI doymamış model → aynı modelde steady akış ile kuplajlı analiz
arasında freatik-yüzey/debi **tutarlı**. Emme işareti: `water_retention` ψ=emme=−(basınç yükü)=y−h
(doygun h≥y → k_rel=1; taban `WaterRetention.k_rel_min`). `retention=nullptr` ise eski **lineer rampa**
(k_min..1, `transition` genişliğinde) birebir korunur (çekirdek testleri regresyonsuz). Doğrulama
(`test_seepage` van-Genuchten dam): Charny debisi **%0.5** (lineerle aynı mertebe), ıslak bölge doygun
kalır (işaret-haritası doğrulandı), Picard 80 adımda yakınsar. **Robustluk dersleri:**
1. **Yakınsama ölçütü YALNIZ doygun bölgede (ψ≥0):** doymamış bölgede k_eff=k·k_min çok küçük → head
   zayıf-koşullu/anlamsız salınır; tüm düğümde ölçmek sahte ıraksama (~6e-3 taban) verir.
2. **Güçlü alt-gevşetme (relax≈0.15):** serbest-yüzey flip-flop'unu (ψ=0'ı geçen düğüm k_rel sıçraması)
   kırar; relax≥0.3 limit-cycle'a takılır, 0.15 ~100 adımda 1e-6'ya yakınsar.
3. **Küçük transition (~eleman boyu) = doğru debi:** geçiş bandı conductance bias'ı debiyi şişirir
   (tr=0.2→%2.1, 0.1→%0.78, **0.05→%0.08**).
**Doğrulama — Dupuit/Charny dikdörtgen baraj** (iki rezervuar h₁>h₂, seepage face yok):
```
q = k(h₁²−h₂²)/(2L)      (CHARNY 1951 teoremi: debi serbest-yüzey şeklinden BAĞIMSIZ, KESİN)
```
*(test_seepage: 73 adım yakınsar, q **%0.08** Charny'ye; freatik yüzey @L/2 Dupuit parabolüne ~%3.6
— Dupuit parabolü yaklaşık, Charny debisi kesin.)* Kaynak: Charny (Polubarinova-Kochina, *Theory of
Groundwater Movement*); Dupuit (1863); Bear (1972).

**Yakınsama ölçütü (kritik):** ham çözüm hn DEĞİL, **GEVŞETİLMİŞ head'in** adımlar arası değişimi
ölçülür: serbest-yüzeydeki tek bir düğüm ham haritada limit-cycle yapar (ψ=0'ı geçerken k_rel sıçraması),
under-relaxation biriktirici head'i yakınsatır ama hn salınmaya devam → |hn−head| sahte tabanda kalır.
Taban ∝ relax × ham-genlik (serbest-yüzey mesh-çözünürlük limiti; ~1e-5 = fiziksel settled).

## 7. Seepage face (sızma yüzeyi) — aktif-set
Serbest-drenajlı sınır (baraj mansap yüzü, kazı yüzü; kuyruk-suyu üstü): su atmosfer basıncında çıkar.
Çıkış noktası a₀ A PRIORI bilinmez. **Aktif-set iterasyonu** (`solve_unconfined_seepage_face`, k_rel
freatik-yüzey iterasyonuyla iç içe; her adımda DofMap yeniden kurulur):
- serbest yüz düğümü ψ=h−y > transition (dead-band) → DEŞARJ (h=y sabitlenir);
- deşarj düğümünde akı içeri (Q>0; outflow Q<0 = gerçek deşarj, Neumann işaret konv.) → serbest.
- **Dead-band (histerezis)** çıkış-noktası yakınındaki belirsiz düğümün (ψ≈0, mesh limiti) exit↔no-flow
  chatter'ını keser. Zor problem → relax≈0.05 (taban tol altına insin).
**Doğrulama — Charny SEEPAGE FACE DAHİL kesin:** dikdörtgen baraj h₁=5, kuyruk-suyu h₂=1, üstte sızma
yüzü → q=k(h₁²−h₂²)/(2L). *(test_seepage: 176 adım yakınsar, q **%0.68** Charny'ye; çıkış noktası
a₀≈1.25 > h₂ doğrulanır.)* Kaynak: Charny; Casagrande (1937) exit-noktası; SEEP/W/PLAXIS review-boundary.

İlgili: [[effective-stress-formulation]] (boşluk→efektif), [[literature-review]], [[analysis-and-structural]].
