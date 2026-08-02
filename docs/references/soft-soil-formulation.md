# Soft Soil Modeli — Kilitli Formülasyon (PLAXIS MMM §10, birebir)

**Birincil kaynak:** PLAXIS 2D 2025.1 Material Models Manual, Bölüm 10 "The Soft Soil model [ADV]"
(metin scratchpad çıkarımından; denklem numaraları manuale). Kullanıcı ilkesi: bu modelin "doğrusu"
tanım gereği PLAXIS'tir — formülasyon birebir alınır. Durum: **AŞAMA 1 (malzeme-noktası çekirdeği)
+ AŞAMA 2 (FE bünye entegrasyonu, §6-7) tamam** (`materials/soft_soil.hpp` +
`material_model.hpp ss_return_core/ss_forward` + `test_soft_soil`); GUI/proje-dosyası + başlangıç
p_p (OCR/POP) + undrained = Aşama 3.

## 1. İzotropik ln-yasa (Denk. 10-5..10-7)
Basınç-pozitif p′ ile (çekirdek içi kabul; çözücü çekme-pozitiftir, dönüşüm HS'deki gibi):
- **Birincil (virgin) yükleme:** ε_v = λ*·ln(p′/p′₀)  (λ* = modifiye sıkışma indisi; ε_v hacimsel
  ŞEKİLDEĞİŞTİRME cinsinden — Burland'ın e-tabanlı λ'sından farkı budur; λ*/κ* oranı λ/κ'ya eşittir).
- **Boşaltma/yeniden yükleme (elastik):** ε_v^e = κ*·ln(p′/p′₀) → teğet hacim modülü
  **K_ur = p′/κ*** (Denk. 10-7; lineer p′-bağımlılık). Elastik sabitler ν_ur + κ* (E_ur/K_ur GİRDİ
  DEĞİL); G = 3K(1−2ν_ur)/(2(1+ν_ur)). p′ birim gerilmenin altına inmez (manual: "minimum value of
  p′ is set equal to a unit stress").
- **Kesin güncelleme:** ln-yasa üstel biçimde tam integre olur: p′_new = p′_old·exp(Δε_v^e/κ*)
  (çekirdek bunu kullanır; sapma-kısmı G(p′) trial-p′'de değerlendirilir — adım-içi donma,
  aşağıda dürüstlük notu).

## 2. Akma yüzeyi: elips cap + MC (Denk. 10-8..10-11)
- **Cap:** f = f̄ − p_p,  **f̄ = q̃²/(M²(p′+c·cotφ)) + p′**, q̃ = σ′₁+(δ−1)σ′₂−δσ′₃,
  δ=(3+sinφ)/(3−sinφ) (HS cap'iyle AYNI q̃, Denk. 10-11 HS §6'ya atıf). f̄=p_p elipsi p′=−c·cotφ
  ile p′=p_p arasında; tepe M-doğrusu üstünde. **Akış cap'te ASOSİYE.**
- **Pekleşme (Denk. 10-10):** p_p = p_p0·exp(Δε_v^p/(λ*−κ*)) (basınç-pozitif ε_v^p ile artar;
  manual: "increases exponentially"). p_p ≥ c·cotφ (c=0'da ≥ birim gerilme) — eşik elipsi.
- **Kırılma:** cap'ten BAĞIMSIZ, standart Mohr-Coulomb (c, φ, ψ; SS varsayılanı ψ=0) — M kritik-durum
  çizgisi DEĞİLDİR (manual açıkça: "failure is not necessarily related to critical state").
- **M ← K0NC (Denk. 10-13, Brinkgreve 1994):** M, birincil 1B sıkışmada K0NC'yi tutturacak şekilde
  seçilir: M = 3·√[ (1−K0NC)²/(1+2K0NC)² + (1−K0NC)(1−2ν_ur)(λ*/κ*−1) /
  ( (1+2K0NC)(1−2ν_ur)·λ*/κ* − (1−K0NC)(1+ν_ur) ) ]. **Bu formülün doğruluğu formülle değil
  DAVRANIŞLA pinlenir:** V&V ödometre birincil yüklemesinde σ′_h/σ′_v → K0NC'yi doğrudan ölçer.

## 3. Parametreler (§10.3)
λ* · κ* · ν_ur (varsayılan 0.15) · c · φ (0 YASAK) · ψ (varsayılan 0) · σ_t · K0NC (→M otomatik).
Alternatif girdi: Cc, Cs, e_init → λ* = Cc/(2.3(1+e)), κ* ≈ 2Cs/(2.3(1+e)) (Tablo 10-2; 2.3 = ln10;
κ* ilişkisi yaklaşıktır — boşaltmada K0=1 varsayımı, manual notu). λ*/κ* tipik 2.5–7.

## 4. Aşama-1 çekirdek şeması ve dürüstlük notları
Asal-uzay, şekildeğiştirme-sürümlü (HS deseninde): (1) üstel-kesin elastik öngörücü; (2) MC kontrol
→ ihlalde mevcut `mc_return_mapping` (SS'in o p′'deki elastik modülleriyle); (3) cap kontrol →
skaler Δλ üzerinde lokal Newton (asosiye akış + üstel pekleşme, D_e dönüşte trial-p′'de DONUK);
(4) köşe: MC-sonrası cap yeniden-kontrol (tek geçiş). NOTLAR: (a) D_e'nin dönüş içinde donması ve
G(p′_tr) seçimi adım-boyu küçüldükçe kaybolan bir integrasyon hatasıdır — V&V bantları buna göre
dürüst (%0.1-0.5 sınıfı; kimlikler [pekleşme yasası, MC kırılma değeri] KESİN test edilir);
(b) teğet Aşama 2'de bağlandı (§7: elastikte analitik, plastikte FD); (c) başlangıç p_p
(OCR/POP, MMM §2.8) Aşama 3'te HS `hs_initial_pp` desenine paralel bağlanır.

## 5. V&V (test_soft_soil; oracle'lar çekirdekle formül paylaşmaz)
(a) izotropik virgin: ε_v(p) = λ*·ln(p/p₀) (200 adım, ≤%0.1) · (b) boşaltma κ* eğimi + p_p sabit ·
(c) yeniden yükleme p_p'yi geçince virgin çizgiye döner (hafıza) · (d) pekleşme yasası KESİN kimlik ·
(e) ödometre birincil: σ′_h/σ′_v → K0NC (M-formülünün davranış piniyle) · (f) drenajlı triaxial:
q_f = MC kapalı-formu (σ₃ sabit mixed-control, HS test deseni).

## 6. Aşama 2a — köşe (kenar) dönüşü: geçerlilik-tabanlı Koiter aktif-küme
q̃ ve MC asal **sıralama** üzerinden tanımlıdır → σ₂=σ₃ (triaksiyel basınç, TC) ve σ₁=σ₂ (TE)
meridyenlerinde akma yüzeyi **koni kenarıdır**. İki ölçülmüş sessiz-yanlış bunu zorunlu kıldı:
1. **Tek-yanlı akış kenarda dallanır:** q̃'nun sıralı ağırlıkları (1, δ−1, −δ) σ₂ ile σ₃'e çok
   farklı yük verir; ödometrik yolda ulp'lik bir yanal asimetri tek-yanlı dala kilitlenip
   **büyür** — iki yanal K0NC'den ayrışıyordu (test g1 yakaladı: σ_xx/σ_v=0.58, σ_zz/σ_v=0.64).
2. **Tolerans-tabanlı rol-ortalama haritayı süreksiz bırakır:** bünye haritası tolerans sınırında
   O(adım) sıçrar → global Newton yakınsayamaz (test g4 yakaladı: residual 2e3'te salınım → NaN).

Tedavi (`return_with_edges`): önce trial sıralamasından **tek-yanlı yüzey dönüşü**; dönüş trial
sıralamasını bozarsa (σ₂<σ₃ vb.) ilgili **kenar dönüşü** — bağlı çiftin rolleri (n_A+n_B)/2 olarak
ortalanır ve dönüş sonrası çift **ortak değere eşitlenir**. w₂=w₃ iken f yalnız σ₂+σ₃'e bağlı
olduğundan eşitleme (saf 2G deviatorik transfer) f'yi ve hacimsel pekleşmeyi değiştirmez.
**İki-çarpanlı kenar çözümüyle özdeşlik:** ayna yüzeyler f_A (roller 1,2,3) ve f_B (2↔3) için
Δλ_A n_A + Δλ_B n_B akışı, Δλ' = Δλ_A+Δλ_B ve t = (Δλ_A−Δλ_B)(n_A2−n_A3)/2 transferiyle tam olarak
"ortalanmış akış + eşitleme"ye eşdeğerdir; yüzey dönüşünün tam kenara indiği sınırda iki harita
birebir örtüşür → **bileşik harita süreklidir, tolerans yoktur** (ölçek/ulp bağımsız).
Sekant Δλ=0 dönerse bu "yüzeyin içindeyiz" değil sekant başarısızlığıdır (tek-yanlı çaprazlanma
kıvrımı f(Δλ)'yı monoton olmaktan çıkarır) → kademe ilerletilir (TC kenarı → tam bağ).
Yüzey gradyanları artık **analitik** (sıralı-rol; cap: a·w+b·1, a=2q̃/(M²p̄), b=(1−q̃²/(M²p̄²))/3)
— Aşama-1'in merkezi-fark gradyanı ve onun gürültü tabanı kalktı.

## 7. Aşama 2b — FE entegrasyonu (`ss_return_core` / `ss_forward`, material_model.hpp)
- **Voigt elastik trial (çerçeve-bağımsız):** üstel-kesin ortalama p_tr = p_c·exp(Δε_v/κ*) +
  trial-p'de G ile lineer deviatorik — izotrop yasa Voigt'ta doğrudan kurulur.
- **Spektral ayrışım + koaksiyel rekonstrüksiyon:** düzlem-içi 2×2 Mohr bloğu + düzlem-dışı normal
  (σ_zz / çember σ_θ) üçüncü asal — HS `hs_return_core` iskeleti birebir (kaynak-rank eşleme).
- **Çekirdeğe TAM TERSİNME ile giriş:** asal Δε = elastik yasanın committed→trial arasındaki tam
  tersi: Δε_v = κ*·ln(p_tr/p_n), Δe = Δs_dev/2G (rank-eşleşmeli koaksiyel — HS'in
  deps_p = C_e(σ_tr−σ_n) hilesinin doğrusal-olmayan yasaya uyarlaması). Bünye TEK kaynaktan
  (`ss_step`) çalışır; elastik adımlar birim dönüşümdür (yuvarlama hatasına dek).
- **Alt-adımlama:** alt-adım başına (|Δε_v| + 2·max|Δε_dev|)/κ* ≤ 0.05 (üstel ortalama ve
  2G ∝ p/κ* aynı 1/κ* ölçeğinde büyür → tek ölçüt). Ölçüm: 0.5 ölçeğinde donuk akış dönüşü
  ıskalıyordu (cap sessizce atlanıp elastik üstel şişme). FD teğetin pertürbe koşuları taban
  koşunun `nsub`'ına SABİTLENİR (HS dersi: serbest ceil sıçraması FD kolonunu kirletir).
- **Teğet:** elastik adımda analitik (K_tr, G) izotrop operatörü (üstel ortalamanın trial'daki tam
  continuum teğeti; dG/dp çapraz terimi donuk-modül geleneğiyle atlanır); plastik adımda İLERİ-FARK
  (h = 1e-6·(1+|ε|); üstel yasanın 1/κ* ölçekli eğriliği yüzünden √eps yerine bilinçli büyük adım).
  SS'te HS'teki "ucuz continuum" yolu YOKTUR: elastik operatör cap üstünde λ*/κ* kat fazla serttir
  → kContinuum/kConsistent her ikisi FD üretir (kNone atlar).
- **Aşama-2 V&V (test_soft_soil g1-g4):** (g1) düzlem-şekildeğiştirme integrate_point ödometre
  zinciri == asal çekirdek zinciri (1e-8) + K0NC tam FE yolundan; (g2) 30° döndürülmüş yol,
  gerilmeler geri döndürülünce özdeş (çerçeve-bağımsızlık, 1e-6); (g3) eksenel-simetrik kol
  (çember=gerçek şekildeğiştirme) aynı çekirdek zinciriyle özdeş (1e-8); (g4) tri15 kapalı kolon
  BVP'si (solve_nonlinear, verilen-oturma, 10 commit'li adım): her Gauss noktası çekirdek
  zincirini ~3e-12 ile üretir, Newton FD teğetle her adımda yakınsar, K0NC BVP içinden gelişir.
