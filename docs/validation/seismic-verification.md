# Sismik/dinamik modül — Doğrulama & Geçerleme (V&V) kaydı

Bu belge KATAI 2D sismik/dinamik izinin (D1–D4a) **her hesabının hangi bağımsız referansa karşı,
hangi toleransla** doğrulandığının savunulabilir kaydıdır. Amaç: bir mühendis/denetçi, KATAI 2D ile
yapılan sismik hesabın doğruluğunu bu kayıttan izleyebilsin. **İlke:** her sonuç kapalı-form veya
BAĞIMSIZ bir yöntemle sınanır; kod tek-kaynaklı öz-tutarlılıkla YETİNMEZ. Sınırlar dürüstçe belirtilir.
İlgili: [[dynamic-seismic-formulation]], [[tbdy-2018-seismic]], [[project-plaxis-parity-roadmap]].

> **Yasal/etik not.** KATAI 2D ile üretilen sismik hesaplar, kullanan mühendisin sorumluluğundadır;
> bu kayıt doğrulama kapsamını şeffaf kılar ama nihai mühendislik yargısının yerini tutmaz. Aşağıdaki
> "Bilinen sınırlar" bölümü henüz doğrulanmamış alanları AÇIKÇA listeler.

## Doğrulanan hesaplar (test_dynamics + test_tbdy_seismic)

| # | Hesap | Bağımsız referans | Sonuç (hata) |
|---|---|---|---|
| D1a | SDOF serbest-titreşim periyodu | T=2π√(m/k) (Chopra) | %0.00 |
| D1a | SDOF sönümlü log-decrement | δ=2πξ/√(1−ξ²) | %0.00 |
| D1a | SDOF zorlanmış kararlı-hal | (F₀/k)·Rd(r,ξ) | %0.03 |
| D1a | Newmark enerji-korunumu (γ=½) | algoritmik-sönüm=0 | genlik korunur |
| D1b | 1D kayma kolonu doğal frekanslar | f_n=(2n−1)Vs/4H (Kramer) | f₁ %0.01 |
| D1c | Taban-harmonik transfer fn (sub-rez.) | 1/cos(ωH/Vs) | %0.1 |
| D1c | Taban-harmonik transfer fn (rezonans) | 2/(πξ) (Kramer Eq 7.30) | %0.0 |
| D2 | 2D tutarlı kütle toplam-korunum | 1ᵀM1=2ρ·alan (tri6/tri15) | ~1e-14 |
| D2 | 2D SH kolon == 1D site-response | f_n=(2n−1)Vs/4H + transfer fn | %0.00 |
| D3 | Absorbing dashpot yansıma katsayısı | E/E0=R²=((Z−c)/(Z+c))² | c=Z→0.0000 |
| D3 | 2D absorbing base radyasyonu | serbest=hapis / absorbing=ışıma | 1.0 / 0.0000 |
| D3b | Free-field yan sınır serbest-alanı korur | u_B==u_A (referans) | 0.0000 |
| D3b | (kontrast) absorbing-only yanlış yutar | 0.9996 (sürücü ŞART) | — |
| D4a | TBDY yatay elastik spektrum kolları | TBDY 2018 §2.3.4.1 (resmî) | KESİN |
| D4a | TBDY zemin katsayıları F_S/F_1 | Tablo 2.1/2.2 (resmî AFAD PDF) | KESİN |
| D4a | Response spectrum (tek sinüs) | A·Rd(r,ξ) kapalı-form | %0.0 |
| D4a | Response spectrum (GENİŞ-BANT) | analitik süperpozisyon (çok-ton) | <%0.4 |
| D4a | Response spectrum (motor-bağı) | inline == solve_newmark | 8.6e-14 |
| D4a | Response spectrum kısa-periyot | sub-stepping vs naive vs A/(2ξ) | naive −19% → sub −4% |
| D3b | 1D free-field kolon (solve_free_field_column) | (4/π)A/(ω₁²·2ξ) rezonans yüzey deplasmanı | %0.1 |
| D3b | (GUI) free-field yan sınır 1D site-response'u geri kazanır | serbest-yan+FF vs 1D teori | −%0.1 (serbest-yan-FF-yok: −%34) |
| D4b | (GUI) 2D dinamik sürücü = 1D site-response (rezonans) | (4/π)A/(ω₁²·2ξ) + peak yüzey ivmesi 2/(πξ)A | %0.6 / %0.2 |
| D5 | Yapısal elastik montaj == statik teğet | doğrulanmış `solve_nonlinear` yolu (BAĞIMSIZ rota) | 9.5e-15 (bağıl) |
| D5 | Rijit yatay öteleme sıfır-enerji modu | K·r=0 (teorem; zemin+plate+interface+ankraj) | 3.8e-17 (bağıl) |
| D5 | Toplam öteleme kütlesi | rᵀMr = ρ·alan + ρA·L (el hesabı) | <1e-10 |
| D5 | Plate konsol öz-frekansları | f_n=(β_nL)²/(2πL²)√(EI/ρA), Blevins Tablo 8-1 | f₁ −%0.01 / f₂ −%0.03 |
| D5 | ρA gerçekten kütle (ölçekleme yasası) | ρA→2ρA ⇒ f→f/√2 | <1e-6 |
| D5 | Gergin geogrid: lineer dal == tension-only dal | statik `solve_nonlinear` | 5.3e-15 (bağıl) |
| D5 | Gevşek geogrid: statik dalda taşımaz | ≡ geogrid'siz model | 1.4e-15 (bağıl) |
| D5 | Kapalı-form başlangıç ivmesi a(0)=−r·a_g(0) | M faktörleyen yol (tekil-olmayan M) | 9.9e-17 |
| D5 | (GUI) yok-olan yapı limiti | EA,EI,w→0 ⇒ zemin-only site-response | +%0.007 |
| D5 | (GUI) gerçek duvar tepkiyi değiştirir | SSI gerçekten kuplajlı (0.4 m diyafram) | +%2.7 |
| D6 | Sismik duvar M zarfı — YARI-STATİK limit | f=f₁/60'ta statik `solve_nonlinear` (−M·r·a_peak) | **+%0.02** |
| D6 | Kuvvet zarfı lineerlik | a_g→2a_g ⇒ kuvvet ×2 (lineer sistem) | KESİN (2.000000000) |
| D6 | (GUI) duvar zarfı üretilir + lineer ölçeklenir | max|M| ×2 | 2.000000 |
| D6b | Arayüz τ zarfı — YARI-STATİK limit | statik `interface_force_diagram` (−M·r·a_peak) | **+%0.04** |
| D6b | Arayüz σ_n zarfı — YARI-STATİK limit | statik `interface_force_diagram` | **+%0.03** |
| D6b | Elastik rapor == Coulomb raporu (elastik dalda) | aynı deplasman alanı, güçlü arayüz | 0.0 (KESİN) |
| D6b | Elastik rapor Coulomb'da KIRPMAZ (zayıf arayüz) | τ_el=4.33 vs τ_Coulomb=0.5 (c_i) | ayrışma KANIT |
| D6b | (GUI) arayüz zarfı üretilir + lineer ölçeklenir | max|τ| ×2 | 2.000000 |
| D7 | **Süperpozisyon KİMLİĞİ (a_g=0)** — duvar momenti | ebeveyn fazın kendi statik çıktısı | **0.000e+00 (KESİN)** |
| D7 | **Süperpozisyon KİMLİĞİ (a_g=0)** — arayüz σ_n | ebeveyn fazın statik σ_n'i (54 kPa, dişi var) | **0.000e+00 (KESİN)** |
| D7 | Talep/kapasite genlikle monotonik | a_g 4.0/1.0/0.2 → %92/%79/%8 aşan | monotonik |

## Yöntemsel güvenceler (kod-hatasına karşı)
- **Bağımsız yöntem çakışması:** response spectrum HEM kapalı-form (tek sinüs) HEM analitik
  süperpozisyon (geniş-bant) HEM `solve_newmark` motoruyla çapraz-doğrulandı. (Bu V&V'nin değeri
  somuttur: geliştirme sırasında yanlış bir Newmark sabiti [a2=2/h yerine 4/h] girildi ve bu üç bağımsız
  kontrol onu ANINDA yakaladı — tek öz-tutarlılık kontrolü kaçırabilirdi.)
- **Sessiz-hata koruması:** `response_spectrum` alt-adımlama yapar (her salınıcı ≥25 adım/periyot);
  kaba-dt'li kayıtta bile kısa-periyot ordinatı SESSİZCE yanlış olmaz (naive tek-adım −19%, sub-stepped
  −4%; yeterli örneklenmiş kayıtta <%1). Kalan hata GİRDİ-örnekleme limitidir (fiziksel), gizlenmez.
- **Enerji tabanlı kontroller:** absorbing/free-field sınırlar enerji ışıması/korunumu ile sınanır
  (yön-bağımsız, işaret-hatasına dayanıklı).
- **Free-field yan sınır (GUI sürücü):** 2D dinamik çözümde serbest (rulo) yan sınırlar sismik taban
  gövde-kuvveti −M·r·a_g altında konsol gibi EĞİLİR → yüzey tepkisi 1D kayma-kolonu site-response'undan
  sapar (ölçülen −%34). Lysmer free-field yan sınırı (kenar dashpot C_b + 1D serbest-alan sürücü kuvveti
  C_b·v_ff; her yan kendi zemin profilinin 1D kolonuyla, free_field.hpp) yan sınırları serbest-alanı
  izlemeye zorlar → 1D site-response'u −%0.1 ile geri kazanır. Doğrulama uçtan-uca, gerçek unstructured
  mesh üzerinde (boundary-edge çıkarımı + 1D kolon + sürücü) `test_dynamic_gui`'de.
- **SSI montajı iki BAĞIMSIZ rotayla çakıştırıldı (D5):** dinamik sistemin yapısal rijitliği, uzun süredir
  doğrulanmış STATİK çözücünün (`solve_nonlinear`; test_plate_soil/test_interface/test_anchor'la pinlenmiş)
  u=0 elastik teğetiyle AYNI olmak zorundadır. Aynı zemin+duvar+arayüz+ankraj problemi iki yoldan çözülür
  ve makine kesinliğinde (9.5e-15) eşleşir → iki DOF-eşleme kuralının sessizce ayrışması imkânsızlaşır.
  Buna ek olarak K·r=0 (teorem) yanlış eşlenmiş bir perde öteleme DOF'unu, rᵀMr=toplam kütle ise plate
  kütlesinin SÜRÜLEN DOF'lara oturduğunu bağımsız olarak sınar. Konsol öz-frekansı (Blevins) ρA'yı gerçek
  bir kütle olarak pinler — rᵀMr tek başına kütle YANLIŞ dağıtılmış olsa da geçerdi (yalnız TOPLAMI sınar).
- **Tekil kütle (SSI'ye özgü sessiz-çökme sınıfı):** yapı dinamik sisteme girince M GENELDE TEKİLDİR
  (plate dönme DOF'ları yalnız ρI taşır; **ağırlıksız plate = GUI varsayılanı w=0** hiç kütle katmaz →
  ölçüldü: 129 denklemin 21'i M'de tam-sıfır satır). Newmark bunu kaldırır (M yalnız matvec; K_eff
  tekil değil) — AMA varsayılan başlangıç ivmesi M'yi FAKTÖRLER. Tekil matrisi PARDISO'nun SPD yoluna
  vermek çözücüyü bozar: erişim ihlali = C++ exception DEĞİL → `/EHsc catch(...)` YAKALAMAZ → uygulama
  mesajsız kapanır (bu projenin iki kez ısırıldığı sınıf). Çözüm yaklaşım değil KESİN: dinlenmeden
  başlarken a(0)=−r·a_g(0), M a(0) = −M·r·a_g(0) = F(0)'ı HER M için (tekil olsa bile) tam sağlar;
  fiziği "t=0'da sistem henüz deforme değil, relatif çerçevede rijit olarak −a_g(0) ile ivmelenir".
  Tekil-olmayan M'de faktörleyen yolla 9.9e-17 örtüşür (yaklaşıklık DEĞİL). `test_ssi_dynamics` (f).
- **Sismik yapısal kuvvet zarfı YARI-STATİK LİMİTLE sınandı (D6):** dinamik run duvarın N/Q/M zarfını
  raporlar (mühendisin kesit tasarladığı sayılar) → dinamik koddan HİÇ pay almayan bir oracle gerekir.
  Sistem ÇOK YAVAŞ sarsılınca (f = f₁/60) atalet+sönüm kuvvetleri elastik olanların yanında yok olur →
  tepki, anlık gövde kuvveti −M·r·a_g(t) altındaki STATİK çözümler dizisine dejenere olur. a_g'nin en
  büyük olduğu anda duvar momenti, −M·r·a_peak altındaki STATİK çözümün (doğrulanmış `solve_nonlinear` +
  `plate_force_diagram`) momentine EŞİT olmalıdır: ölçülen **+%0.02** (kalan artık atalet ~(f/f₁)²≈3e-4).
  Bu tek kontrol tüm zinciri pinler: dinamik deplasmanlar, onlardan kuvvet geri-kazanımı, ve zarf birikimi
  (zarf sessizce peak-DEPLASMAN anını raporlasa ya da bir istasyon düşse bu sayı kayar). `test_ssi_dynamics` (g).
- **ARAYÜZ SİSMİK RAPORU ELASTİK OLMAK ZORUNDA (D6b — sessiz tutarsızlık yakalandı):** statik
  post-processor `interface_force_diagram` Coulomb return UYGULAR (τ'yu dayanımda kırpar, tension cut-off,
  σ_n0 ekler). Dinamik dal ise arayüzü **elastik (k_n,k_s)** çözer ve sıfırdan başlar. O post-processor'ı
  dinamik fazda OLDUĞU GİBİ kullanmak, çözücünün ÜRETMEDİĞİ bir gerilmeyi raporlardı = kendi
  deplasmanlarıyla DENGEDE OLMAYAN çıktı (sessiz tutarsızlık). ÇÖZÜM: `interface_force_diagram(..., elastic)`
  → τ = k_s·Δu_s, σ_n = k_n·Δu_n (kapak yok, cut-off yok, σ_n0 yok = dinamik sistemin gerçekten çözdüğü
  bünye). Test (h) bayrağın KENDİSİNİ iki yönden pinler: güçlü arayüzde (Coulomb dalı ısırmaz) iki rapor
  BİREBİR aynı (0.0 → bayrak sapma katmıyor), zayıf arayüzde (c_i=0.5) ayrışma GERÇEK ve ölçülü.
- **✅ D7 — TOPLAM TASARIM ETKİSİ + TALEP/KAPASİTE (aşağıdaki "ölçülen risk" ARTIK ÖLÇÜLÜYOR):** dinamik faz
  artık ebeveyn fazın STATİK durumunu süperpoze edip **TOPLAM** (statik+dinamik) tasarım etkisini raporluyor
  (PLAXIS de toplamı raporlar, Ref §9.4.5). **İŞARETLİ min/max zarf ŞART oldu:** statik ofset simetriyi
  bozar → max_t|N_s+N_d(t)| = max(|N_s+min N_d|, |N_s+max N_d|); max|N_d| yalnız N_s=0'da doğrudur.
  **DOĞRULAMA — SÜPERPOZİSYON KİMLİĞİ:** a_g=0'da dinamik artım özdeş sıfırdır → raporlanan toplam,
  ebeveyn fazın statik kuvvetlerine **BİREBİR** eşit olmalı: ölçülen **0.000e+00** (hem duvar M'si hem
  arayüz σ_n'i; σ_n statik 54 kPa → kontrolün dişi var). Bu kimlik, istasyon eşlemesi/işaret/düşen terim
  hatalarının hiçbirini sağ bırakmaz. Süperpozisyon LİNEER artım için KESİNDİR: K u_s = f_s ve
  M ü_d+C u̇_d+K u_d = −M·r·a_g toplanır (u_s zamandan bağımsız). **TALEP/KAPASİTE:** τ_max toplam σ_n'e
  bağlı olduğundan ancak statik durum gelince hesaplanabilir → `utilisation = |τ_total|/τ_max` istasyon
  başına + `max_utilisation` + `over_fraction`. **VERDICT DEĞİL ÖLÇÜ:** yüzeyde σ_n→0 olduğundan τ_max→c_i
  ve HER duvarın tepesinde kısa bir bölge gerçek sarsıntıda kapasiteyi aşar — bu gerçek yerel kaymadır,
  tüm analizi mahkûm etme gerekçesi değil → hem tepe oran hem AŞAN YÜZDE veriliyor, karar mühendisin.
  τ_max=0 (arayüz ayrılmış = hiç kayma kapasitesi yok) → oran sınırsız → `kUtilCap=999` ile kapanır.
  **DERS (test öncülüm çürüdü): `Rinter` ile sıralama YAPILAMAZ** — PLAXIS'in sanal-kalınlık formülasyonu
  rijitliği Rinter² (G_i=Rinter²G → k_s=G_i/t_i), dayanımı Rinter¹ düşürür → zayıf arayüz aynı zamanda
  YUMUŞAKTIR, τ=k_s·Δu_s kapasiteden hızlı düşer, utilisation ~ Rinter → **zayıf arayüz DAHA DÜŞÜK oran
  verir** (ölçüldü: Rinter 0.05 → 0.99× vs Rinter 1.0 → 3.42×). Ayrımcılık GENLİKLE sınanır (monotonik).
  **DERS-2: `rinter_rigid` varsayılanı TRUE ve `Rinter`'i EZER** (`Rinter = rigid ? 1.0 : R`) — testte
  temizlenmezse arayüz rijit kalır (bu, ilk testimin sessizce aynı sonucu vermesinin nedeniydi).
- **🔎 (TARİHÇE) ÖLÇÜLEN RİSK — arayüz kayma talebi (D7 ile kapatıldı):** test (h)'de zayıf arayüzde elastik
  talep **τ=4.33 kPa**, gerçek Coulomb kapasitesi **0.5 kPa** → lineer dinamik analiz kapasitenin ~**8.7
  katını** taşıtıyor ve sonuç orada GÜVENLİ TARAFTA DEĞİL. Bu, lineerleştirmenin en tehlikeli yüzü.
  ŞU AN: sınır τ'nun okunduğu HER yere yazılıyor (Output→Interfaces sarı uyarı+help, rapor 6b notu EN/TR,
  faz sonuç mesajı) ve `[bonded]` rozeti dinamik fazda GÖSTERİLMİYOR (`any_slip=false` arayüz yapışık
  olduğu için değil, KAYAMADIĞI için → "[elastic, no slip check]").
  **SONRAKİ (ayrı, dikkatli iş): talep/kapasite oranı.** Rijit gerekçe: dinamik tepki sıfır etrafında
  simetrik olduğundan max_t|τ_statik+τ_dyn(t)| ≥ max_t|τ_dyn(t)| (bir anda işaretler çakışır) →
  **|τ_dyn|_max > τ_max ⟹ gerçek arayüz KESİNLİKLE kayar ⟹ lineer sonuç orada geçersiz** (TEK YÖNLÜ
  kanıt; oran < 1 geçerlilik KANITLAMAZ çünkü statik kayma eklenir). Tam/rigorous kontrol τ_max'ın toplam
  σ_n'e bağlı olmasını ve statik durumun dinamik faza TAŞINMASINI gerektirir (aşağıdaki "statik durumdan
  başlamaz" sınırı) → yarı-proxy uydurmak yerine ayrı artım olarak bırakıldı.
- **Zarf ≠ denge durumu (raporlamada dürüstlük):** dinamik faz kuvvetleri max|·| ZARFIDIR — farklı
  istasyonların uç değerleri genelde FARKLI anlarda oluşur, dolayısıyla zarf bir denge durumu değildir.
  Ayrıca yalnız DİNAMİK etkidir (faz dinlenmeden başlar) → toplam tasarım etkisi için statik faz kuvvetleriyle
  süperpoze edilmelidir. Zarf olduğu `StructForce::envelope` ile taşınır ve **hem GUI'de (sarı uyarı + help)
  hem HTML raporunda (bölüm 6 notu, EN/TR)** yazılır; `.res` formatı v3'e yükseltildi (v1/v2 dosyalar
  sismik-öncesi → envelope=false doğru). Tek anlık snapshot raporlamak duvar momentini EKSİK gösterirdi
  (peak moment, peak deplasman anında oluşmaz).

## 🔍 SİSMİK DENETİM (2026-07-16) — bulunan ve DÜZELTİLEN sessiz-yanlış kusurlar
Kullanıcı direktifi ("sunulan her sonucun doğru olduğundan emin olmadan ilerlemeyiz") üzerine modül
uçtan uca denetlendi + PLAXIS 2D 2025.1 resmî manualleriyle (Ref/Sci/Mat/Tut) karşılaştırıldı. Bulunan
kusurların HEPSİ "makul görünen yanlış sayı" sınıfındaydı — çökme değil, sessiz hata. Hepsi düzeltildi
ve her biri için REGRESYON TESTİ yazıldı (`test_dynamic_gui` A1-A6).

| # | Kusur | Etki (ölçülen) | Düzeltme |
|---|---|---|---|
| A1 | **Dinamik faz staged-construction maskesini YOK SAYIYORDU** — `assemble_stiffness`/`assemble_mass`'te `active_element` parametresi yoktu (kardeş `assemble_gravity`'de vardı) | Kazılmış zemin K ve M'ye tam katkı veriyordu; GUI kazıyı ÇİZİYOR, çözücü dolu modeli çözüyordu → duvar sismik kuvvetleri EKSİK = **güvensiz tarafta**. Ölçülen: 10 m kolondan 4 m kazınca yüzey ivmesi **2.45 → 9.03 m/s² (%269)**; düzeltmeden önce bu iki koşu BİREBİR AYNIYDI (f₁ 1.98→3.30 Hz kayması tamamen kayboluyordu) | `active_element` çekirdek API'ye eklendi (varsayılan `{}` = hepsi aktif → geriye uyumlu); dinamik dal `act` geçiriyor |
| A2 | **HS/HSsmall'da dinamik faz GİRİLMEMİŞ E ile çalışıyordu** — malzeme editörü E'yi yalnız LE/MC'de gösteriyor (`main.cpp:3408`), HS kullanıcısı E50/Eoed/Eur giriyor → `M.E` struct varsayılanı 1.3e4'te kalıyor | TÜM sismik koşu keyfi 13 MPa üzerinden: Vs, f₁, büyütme, duvar kuvvetleri modelle İLGİSİZ. Sessiz | Dinamik fazda HS/HSsmall açık mesajla REDDEDİLİYOR (PLAXIS Mat.Man §3.5 aynı şeyi kullanıcıdan ister: dinamikte rijitlik doğru dalga hızını vermeli = küçük-şekildeğiştirme rijitliği) |
| A3 | **`E'_inc` / `c'_inc` (derinlikle artan rijitlik/dayanım) editörde var, yardım metniyle vaat ediliyor, dosyaya kaydediliyor — ÇEKİRDEKTE HİÇ KULLANILMIYOR** (`E_inc` yalnız main.cpp/project.hpp/project_io.hpp'de) | **HER FAZ** sessizce üniform E'_ref ile çözüyordu. Sismikte E(y) = Vs profilidir → en yaygın gerçek girdi yok sayılıyordu. En GENİŞ kusur (dinamiğe özgü değil) | `unsupported_profile_warning` → TÜM fazlarda erken red (gradyan verilmişse). Gerçek E(y) implementasyonu ayrı iş |
| A4 | Dinamik kütle her yerde `gamma_unsat`; `gamma_unsat=0` → `assemble_mass` her elemanı atlıyor → F=−M·r·a_g=0 | Su altında atalet **%15 eksik** (güvensiz); sıfır ağırlıkta `ok=true, peak\|u\|=0` = **makul görünen sıfır** | Yeni `assemble_mass_phreatic` (Gauss-noktası başına: su altında ρ_sat, üstünde ρ_unsat — `assemble_gravity_phreatic` ile AYNI desen); sıfır kütle (rᵀMr=0) erken red |
| A5 | Yüzey monitörü tüm mesh düğümlerini tarıyordu (aktiflik/serbestlik filtresi yok) → kazılmış modelde pasif bloğun orphan (pinlenmiş) düğümüne düşebiliyordu; `surf_eq<0` → `a_surf = a_g(t)` | Raporlanan "yüzey response spectrum" **GİRDİ spektrumunun kendisi** oluyordu — TBDY tasarım spektrumunun üstüne saha tepkisi diye çiziliyordu | Monitör yalnız AKTİF elemana ait + ux'i SERBEST düğümler arasından seçiliyor; hiçbiri yoksa erken red (asla a_g yayınlanmıyor). Test: kazılmış model 18.49 m/s² (girdi PGA 1.0 DEĞİL) |
| A6 | Taban hareketi `a_g(k·dt)` ile noktasal örnekleniyor, anti-aliasing yok; Newmark koşulsuz kararlı → az adım sessizce BAŞKA bir sinyali entegre ediyor | **GUI VARSAYILANI** (1.0 s / 25 adım @ 2 Hz) = 12.5 adım/çevrim — programın KENDİ yardım metni "20-40" diyor. 10 s/50 adım @ 5 Hz = 1 adım/çevrim = Nyquist altı, hatasız yanlış cevap | `dynamic_step_warning` (konsolidasyon Δt_crit deseniyle aynı) → faz editöründe sarı uyarı; Nyquist altını açıkça söylüyor; Ricker ~2.5f üzerinden değerlendiriliyor |
| A7 | Free-field yan sınır seçimi `\|nx\|≥0.5` + `extract_boundary_edges` (tek aktif elemana ait HER kenar) | Duvar DİKİŞİ (split ikiz düğüm → iki yüz de sınır görünüyor), KAZI yüzü ve 30°'den dik HER ŞEV yan sınır sanılıyordu → model ORTASINA dashpot+sürücü; sol/sağ gruplaması gerçek yan zinciriyle iç zinciri TEK kolona birleştiriyor (örtüşen y + yüzen ikinci dal) → çöp v_ff, hatasız | Kenarın TÜM düğümleri modelin gerçek sol/sağ uç düzleminde olmalı (`on_extreme`) → yalnız hakiki uzak-alan yanları |

**PLAXIS KARŞILAŞTIRMASI (resmî manuallerden, denetimle birlikte yapıldı) — mimari sapma kayda geçti:**
PLAXIS'te dinamik faz **ebeveyn fazın durumundan DEVAM EDER** (Sci §6.4: sınır gerilmesi başlangıç hızına
çevrilir, "önceki hesap veya başlangıç gerilme durumuna dayanır"; Ref §7.9.2.3: "*Reset displacements to
zero* **gerilme alanını ETKİLEMEZ**" = yalnız deplasman datumu) ve **TAM NONLİNEER** çalışır (Ref §11.10.4:
Dynamics, plastik analizle AYNI yakınsama kriterlerini kullanır; Sci §6.1: "prensipte PLAXIS'teki tüm
modeller dinamik analizde kullanılabilir"). Bizim dinamik fazımız **sıfırdan başlayan LİNEER** bir
analizdir → bu, A1/A2'nin de kökü olan bilinçli mimari sapmadır ve aşağıda "bilinen sınırlar"da kayıtlı.
PLAXIS zarfı **İŞARETLİ min/max** çiftidir (Ref §9.4.5: "historical maximum and minimum forces"), bizimki
max|·| — lineer + dinlenmeden başlayan sistemde tepki sıfır etrafında simetrik olduğundan max|·| tasarım
değeridir, ama statik durum taşınırsa İŞARETLİ zarf ŞART olacak. PLAXIS raporladığı kuvvet **TOPLAM**
(statik+dinamik); bizimki yalnız dinamik → süperpozisyon kullanıcıya bırakılıyor (yazılı).
**⚠️ KENDİ İDDİAMIN DÜZELTİLMESİ:** D5'te plate dönme ataleti ρI=ρA·d²/12 "PLAXIS Ref §5.6 eşdeğer kalınlık"
gerekçesiyle sunulmuştu. Araştırma: **hiçbir PLAXIS manuali (Ref/Sci/Mat/Tut) plate/kiriş mass matrisinde
dönme ataletinden söz ETMİYOR** — Sci'nin eleman bölümü yalnız RİJİTLİK matrisini veriyor; PLAXIS varsayılanı
**lumped** kütledir (Ref §7.9.5.2, 0..1 arası ayarlanabilir) ve standart FE pratiğinde lumped Mindlin kirişi
yalnız ÖTELEME kütlesi taşır. Yani ρI **bizim mühendislik kararımızdır, PLAXIS paritesi DEĞİLDİR** (d = eşdeğer
kalınlık tanımı PLAXIS'ten, ρI'nin kendisi değil). Konsol öz-frekansı doğrulaması (Blevins, %0.01) geçerli
kalır — o kapalı-form teoriye karşıdır, PLAXIS'e karşı değil. PLAXIS'in `w` [kN/m/m] → ρA=w/g yorumu ise
manuallerde tek cümlede yazmıyor; araştırma bunu Tut §7.6.1 diyafram duvarı sayılarıyla (w=8.3, d_eq=0.346 →
γ≈24 kN/m³) sayısal olarak doğruladı.

## 🔗 TRACK 1a (2026-07-19) — ebeveyn YAPISAL durumunun nonlineer dinamik artıma taşınması (2e42112, eaacb4c)

Nonlineer dinamik çözücü zemin gerilmesini ebeveynden devralıyordu ama YAPISAL elemanlar (arayüz
kayması, ankraj U_p, geogrid ε_p, gömülü-kiriş skin/foot) SIFIR statik ön-yük durumundan başlıyordu →
Coulomb/akma/gevşeme kapakları yalnız dinamik artıma uygulanıyor, statik ön-yük olan her yerde
kayma/akma BAŞLANGICI yanlıştı (D7'nin 8.7× ölçtüğü sınıf). Artık `StructuralInit` ile ebeveynin
yakınsamış deplasman DATUMU + committed plastik durumu taşınır; kapaklar TOPLAM (statik+dinamik)
etki üzerinde çalışır. V&V (bağımsız oracle):
- **Kapasite CEBİRİ (test_dynamics_nonlinear (d)):** 0.6·F_max'a ön-yüklü ankraj, ±0.6·F_max
  yarı-statik sarsıntıda YALNIZ taşıma varken akar (toplam 1.2·F_max; artım tek başına elastik,
  U_p bit-birebir 0); toplam 0.85·F_max kalınca taşımayla bile AKMAZ (U_p bit-birebir == tohum) →
  başlangıcı taşıma mekanizması değil, TOPLAMIN kapasiteyi kesmesi yönetiyor.
- **Sıfır-kuvvet kimliği:** tohumlu + sürücüsüz çözücü TAM hareketsiz kalır (peak |u| tam 0) ve
  committed durum bit-birebir tohuma eşit (datum + tohum + f_int₀ baseline karşılıklı tutarlı).
- **Uçtan uca kimlik (test_dynamic_gui):** K0 → yüklü Plastic → sürücüsüz nonlineer Dynamic
  zincirinde her arayüz VE duvar istasyonu ebeveynin statik değerine **0.000e+00** ile eşit
  (diş: σ_n=54 kPa, M=6.03 kNm/m) — D7 süperpozisyon kimliğinin çözücü-İÇİ hali.
- **Utilisation geri geldi (tutarlı):** taşıma varken talep/kapasite ANLIK (aynı adımın τ ve σ_n'i)
  birikir → çözücü kapağı toplamda uyguladığı için oran yapısı gereği ≤ ~1; over_fraction ≈ 0
  (lineer yolun yalnız ÖLÇTÜĞÜ aşım artık ÇÖZÜLÜYOR). Taşınamayan ebeveynde (boyut uyuşmazlığı /
  diskten yüklenmiş sonuç / Dynamic ebeveyn) dürüst geri-düşüş + faz mesajında açık beyan.
- **any_slip / yielded artık KESİN durum-olayı:** committed plastik durumun tohuma göre bit-düzeyi
  değişimi (elastik dal committed değeri bit-birebir döndürür). Eski `|du_s|>1e-12` türetimi her
  YAPIŞIK eklemi [SLIPPING] damgalıyordu (elastik eklemde du_s=τ/k_s her zaman ≠0) — düzeltildi;
  commit-sonrası bünye yeniden-değerlendirmesi tam akma sınırında eşitlik yazı-turasıydı — kullanılmıyor.
- **Rozet dürüstlüğü (.res v5, `slip_checked`):** lineer zarf "[elastic, no slip check]" der;
  nonlineer zarf gerçek [SLIPPING]/[bonded] bulgusunu taşır (3 UI yüzeyi: panel, HTML rapor, metin).
- **Lineer zarf ankraj/geogrid tutarlılığı (eaacb4c, D6b kuralının yapısal karşılığı):** lineer
  sistem ankrajı KAPAKSIZ, geogrid'i tam-EA çözer → zarf da öyle raporlamalı. Ölçüldü: F_max =
  talep/2 iken lineer zarf 156.5 kN/m (önceden sessizce 78.2 = F_max) + akma iddiası yok; nonlineer
  zarf tam F_max'ta + yielded=yes. Kapaklı diyagram talebi EKSİK raporlar = güvensiz taraf.

## 📊 BENCHMARK Dalga-1 (2026-07-19) — çok-tabaka site-response ✅ (Track 7 ilk kalemi)

**KATAI 2D dinamik GUI yolu, SHAKE-sınıfı KATMANLI viskoelastik transfer-matris çözümünü üretiyor**
(`test_site_response_benchmark` + `docs/validation/site-response-benchmark.md`). Oracle bağımsız
(kompleks transfer-matris; çözücüyle sıfır kod paylaşımı) ve FE'nin Rayleigh sönümünün TAM sürekli
karşılığıyla (G*=G(1+iωβ), ρ*=ρ(1−iα/ω)) kurulu → like-for-like. Oracle'ın kendisi iki kapalı-formla
pinli (1/cos(ωH/Vs) 6e-16; iki-tabaka karakteristik denklemi tan·tan·α=1.0005). Sonuç: tek-tabaka
%0.05/−0.18/%1.7; **iki-tabaka (empedans kontrastı 0.343): flank +%0.06, 1. rezonans tepe −%0.09
(|T|=17), çukur +%0.65, 2. rezonans +%0.38.** Fizik dersi: kuvvetli kontrastta çeyrek-dalga
seyahat-süresi kuralı f₁'i %22 alçak tahmin eder (2.34 vs gerçek 3.005 Hz) — Rayleigh hedef
frekansı seçiminde mühendise not.

## Bilinen sınırlar (henüz DOĞRULANMADI — dürüst beyan)
Aşağıdakiler v1.0 kredibilitesi için tamamlanmalı; şu an KAPSAM DIŞI:
- **PLAXIS 2D Dynamics kıyası YOK.** Çekirdek kapalı-form doğrulandı ama PLAXIS ile nicel parite henüz
  koşulmadı (roadmap ilkesi: PLAXIS parite). 
- **SHAKE/DEEPSOIL kıyası ◑:** LİNEER çok-tabaka büyütme ✅ (yukarıdaki Dalga-1; SHAKE'in lineer
  çekirdeği aynı transfer-matris yöntemidir). KALAN: eşdeğer-lineer iterasyon (G/G₀−γ, ξ−γ eğrileri),
  gerçek kayıt girişi, outcrop/deconvolution (absorbing tabana bağlı), frekans-bağımsız histeretik
  sönümle kıyas.
- ✅ **Gerçek deprem kaydı VAR (akselerogram girişi + El Centro 1940 NS; real-record-elcentro.md):**
  kayıt kimliği yayımlanmış PGA'ya birebir (0.31882 vs 0.319 g), %5 spektrum yayımlanmış büyütme
  bandında (2.87×PGA), gerçek-kayıt ürün koşusu rijit-vs-compliant fiziksel yönde (1.59 vs 0.52 g).
  AÇIK kalan: nokta-nokta yayımlanmış spektrum TABLOSU kıyası (bantlar dürüst; motor kapalı-formla pinli).
- **AFAD Türkiye Deprem Tehlike Haritası web servisi ile uçtan-uca kıyas YOK** (S_S/S_1 girdisi → spektrum;
  formül+tablolar resmî PDF'ten doğrulandı ama servis çıktısıyla birebir kıyas yapılmadı).
- **Kapsam:** yalnız YATAY SH (uy=0). Düşey (P-SV) spektrum §2.3.5, çok-tabaka, elastoplastik dinamik
  (histeretik sönüm), sıvılaşma → sonraki fazlar. Rayleigh sönüm frekans-bağımlı (bant-dışı fazla söner).
- **SSI lineerleştirmesi artık yalnız VARSAYILAN (lineer) dalda.** Lineer dinamik sistemde yapısal
  elemanlar gergin/kaymayan durum etrafında lineerleştirilir: geogrid basınçta da EA taşır, ankraj
  akmaz, arayüz kaymaz (en keskin durum gevşek geogrid, ölçülen %63 — `test_ssi_dynamics` (e));
  hangi basitleştirmenin uygulandığı faz mesajında YAZILIR ve zarf raporu da aynı elastik bünyeyi
  kullanır (eaacb4c: kapaksız ankraj + tam-EA geogrid — kapaklı rapor talebi eksik gösterirdi).
  **Opt-in NONLİNEER dal (Phase.dynamic_nonlinear) bunların hepsini ÇÖZER:** zemin plastikleşir,
  arayüz kayar, ankraj akar, geogrid gevşer — ebeveyn yapısal durumu üzerinden (Track 1a).
- **🏛️ (ESKİ) MİMARİ SAPMA — BÜYÜK ÖLÇÜDE KAPANDI.** PLAXIS gibi (Sci §6.4; Ref §7.9.2.3, §11.10.4):
  nonlineer dinamik faz ebeveynin GERİLMESİNDEN (init_states, c00fa1d) ve YAPISAL DURUMUNDAN
  (Track 1a, 2e42112) devam eder, tam nonlineer çözülür; işaretli min/max zarf + statik süperpozisyon
  D7'den beri var. KALAN: varsayılan LİNEER dal tasarım gereği artım-tabanlıdır (süperpozisyon D7
  post ile) — bit-birebir hız yolu; ve absorbing taban aşağıda ayrı madde.
- ✅ **Gömülü kiriş dinamikte (cc96be1).** Mesh-uyumsuz skin kuplajı eleman-tipine dispatch ile
  çözüldü; K·r=0 kimliği (3.83e-17) ve rᵀMr kütle kimliği kanıt. (Eski erken-red kalktı.)
- ✅ **Arayüz sismik kayma: nonlineer dal ÇÖZÜYOR (4faf467 + Track 1a).** Lineer dal talep/kapasite
  ORANINI ölçer ve yazar (8.7× sınıfı uyarı); nonlineer dal Coulomb kapağını TOPLAM etki üzerinde
  uygular → gerçek kayma raporlanır, utilisation ≤ ~1, over_fraction ≈ 0.
- ✅ **Sismik zemin gerilmesi nonlineer dalda geri-kazanılıyor (7cbae24).** Deprem-sonrası committed
  σ alanı `recover_nodal_stresses_from_gauss` ile raporlanır. Lineer dal sıfır-gerilme datumunda
  çözdüğü için orada ertelendi (bilerek).
- **Plate ağırlığı varsayılan 0'dır (w=0).** Ağırlığı girilmemiş bir duvar dinamikte KÜTLESİZ bir
  sertleştirici gibi davranır (rijitlik katkısı var, ataleti yok) — fiziksel olarak eksiktir. GUI faz
  editörü bunu yazar; kullanıcı plate malzemesinde w'yi girmelidir.
- ✅ **COMPLIANT (absorbing) TABAN VAR (opt-in `Phase.seismic_compliant_base`; formülasyon
  dynamic-seismic-formulation.md §11'de birebir kaynak alıntılarıyla kilitli).** Joyner-Chen: taban
  Lysmer ρV_s dashpotları (en derin tabaka = yarı-uzay onu sürdürür) + girdi AYNI dashpot matrisi
  üzerinden 2ρV_s·v_up traksiyonu (v_up = within hareketinin YARISI, Tut §17.8.5); çözüm TOPLAM
  harekette. **V&V (test_compliant_base, radyasyon-BC'li bağımsız oracle):** oracle kapalı-forma
  2e-16; üniform kolon (mükemmel radyasyon, |T|~1) FE −%0.02; iki-tabaka tepe −%0.05; **DİŞ: rijit
  taban aynı profilde radyasyon çözümünün 8.7 KATI çınlıyor** (kaldırılan kutu-rezonans hatası).
  **Kombinasyonlar da doğrulu:** nonlineer×compliant lineer-limit kimliği **6e-10** (toplam-hareket
  baseline argümanı nonlineer çözücüde birebir tutuyor) + zayıf Tresca kolonu akışa geçiyor (%1.5);
  free-field yanlar×compliant (1D yan kolonlar da compliant tabanlı — `solve_free_field_column_
  compliant`) 1D radyasyon kimliğini **−%4.8** ile korur (SH yanlarda −%0.02; fark köşe-dashpot
  örtüşmesi + serbest-yan 2D etkisi — PLAXIS'te de aynı köşe yapısı; bant %5, dürüstçe kayıtlı).
  Kapsam: yatay SH (taban u_y sabit; mesajda beyan). Nonlineer Δt uyarısı: nonlineerde doğruluk+
  yakınsama dili eklendi (uydurma ikinci eşik YOK; ölçülü kriter refinement çalışması olarak kayıtlı).

## Reprodüksiyon
`cmake --build build/msvc-rwdi --target test_dynamics test_tbdy_seismic test_ssi_dynamics test_dynamic_gui test_dynamics_nonlinear test_site_response_benchmark`
→ çalıştır. SSI (D5) çekirdek montajı: `test_ssi_dynamics` (MKL'siz, pure-Eigen); GUI-yol SSI: `test_dynamic_gui`.
Formülasyon: `dynamic-seismic-formulation.md` (D1–D3b), `tbdy-2018-seismic.md` (D4a, resmî TBDY
değerleri). GUI-yol sürücü + free-field yan sınır: `test_dynamic_gui` (build_problem Dynamic dalı).
