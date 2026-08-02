# PLAXIS 2D ↔ KATAI 2D — Kapsamlı Karşılaştırma Tablosu

**Tarih: 2026-07-20.** Bu tablo, KATAI 2D'nin PLAXIS 2D'ye göre nerede durduğunun SAYISAL kaydıdır.
Dört kaynak sınıfı: (A) **resmî PLAXIS Validation Manual vakaları birebir koşuldu** (üç kolon:
analitik | PLAXIS'in yayımladığı | KATAI); (B) PLAXIS'in doğrulama setiyle **aynı analitiğe** iki
tarafın bağımsız doğrulaması; (C) **doğrudan PLAXIS-sonucu** kıyasları (tutorial + hakemli makale);
(D) dinamik/sismik (PLAXIS bu sınıfta sayı yayımlamaz → bağımsız oracle + yöntem paritesi).
Her satır bir ctest regresyonudur ya da docs/validation kaydına bağlıdır. Dürüstlük notları en altta.

## A. Resmî PLAXIS Validation Manual vakaları — birebir koşu (`test_plaxis_validation`)

Kaynak: PLAXIS Validation Manual (V8; şirketin kendi yayımladığı benchmark seti, analitik + kendi
sonucu birlikte verilir).

| Vaka | Analitik | PLAXIS (yayımlı) | KATAI 2D | KATAI vs analitik | KATAI vs PLAXIS |
|---|---|---|---|---|---|
| §2.1 Pürüzsüz rijit şerit temel, elastik (Giroud) — F [kN/m] | 15.15 | 15.24 (+%0.6) | **15.35** | +%1.4 | **+%0.8** |
| §2.2 Gibson zemini E=299z, ν=0.495 — oturma [m] | 0.0500 (yarı-uzay) | 0.0470 (−%6, sonlu tabaka) | **0.0451** | −%9.7 | **−%4.0** (aynı sonlu-tabaka sapması) |
| §3.1 Dairesel temel taşıma, MC axisym (Cox) — p_max [kPa] | 225.6 | 220.0 (−%2.5) | **234.4** | **+%3.9** | +%6.5 |
| §3.2 c(z) kilinde şerit temel (Davis-Booker, pürüzsüz) — p_max [kPa] | 7.80 | 7.86 (+%0.8) | **8.02** | **+%2.8** | +%2.0 |

Koşu notları (kayıtlı dersler): §3.1'de Cox slip-line çözümü **asosiye** limit yüktür → ψ=φ ile
koşuldu (ψ=0 kıyası bir modelleme farkını doğrulama sayısına karıştırırdı); 0.5 m mesh çökmeyi +%9
fazla tahmin etti (temel yarıçapında 2 eleman — PLAXIS'in kendi düşük-derece uyarısının tri15'teki
hafif hali), 0.25 m ile +%3.9. §2.2'de yüzeyde E→0 uyumu ince yüzey bandında toplar; 0.5 m mesh
−%7, 0.125 m ile PLAXIS'e −%4 (ν=0.495 hafif hacimsel kilitleme payı dahil).

## B. Aynı analitiğe iki taraflı doğrulama (problem sınıfı eşlemesi)

| Doğrulama sınıfı | PLAXIS'in yayımlı sonucu | KATAI'nin sonucu (test) |
|---|---|---|
| Şerit temel taşıma, pürüzlü c(z) (D&B §3.2) | 9.25 vs 9.1 (+%1.6) | (pürüzsüz varyant A'da +%2.8; pürüzlü koşu kuyruk) |
| Prandtl/Reissner taşıma gücü | (V&V setinin temel sınıfı) | Nc=5.200 vs 5.142 (**+%1.1**, tri15 120 elem); φ=20°: 15.13 vs 14.84 (+%2.0) |
| Kiriş eğilmesi (§2.3): u_max | 13.96 vs 13.89 mm (+%0.5) | Konsol kapalı-formu **%0.01** (Blevins); kuartik plate 2.5e-13 |
| Silindirik kavite genişlemesi (§3.4) | eğri üstünde ("very well") | Tam-plastik silindir p=2c·ln(b/a) **%0.3**; Lamé elastik 1e-5 |
| Arayüz Coulomb kayması (§3.3 blok) | 60.4 vs 60 kN/m (+%0.7) | Coulomb return kapalı-form testleri (tam; test_interface) + K0-seed denge kimlikleri |
| 1B konsolidasyon (§4.1 Terzaghi) | "analitiğe yakın" (eğri) | U-T_v eğrisi ≤%2 (test_consolidation/_gui) |
| Serbest-yüzeyli sızma (§4.2 Dupuit) | 0.152 vs 0.150 (+%1.3) | Dupuit + Charny sızma-yüzü + Darcy kapalı-formları (test_seepage) |
| Perde etrafında basınçlı akış (§4.3 Harr) | 0.818 vs 0.80 (+%2.3) | Baraj-altı akış-ağı/uplift kapalı-formları (test_seepage) |
| Şev güvenliği (φ-c reduction) | (V&V pratiği Bishop kıyası) | FoS 1.010 vs Bishop ~0.99 / T6 0.997 (**+%1.3–2.0**) |
| MC çekme kesmesi (MMM Denk. 3-11, §3.3.10 hendek) | (PLAXIS varsayılan açık, σ_t=0) | 6 bölge tersinmesi **kesin** + 8000'lik uygunluk taraması; hendek BVP kret σ₁: kapak AÇIK 0.40 / KAPALI 7.71 kPa (test_mohr_coulomb, test_input_audit) |
| Plate M-N plastik mafsalı (MMM §18.3 elması) | (elmas diyagram + gerilme-noktası kontrolü, Fig 18-1) | Dönüş haritası **CPP-kesin** (81×81 tarama f=0.00, brute-force enerji-normu oracle birebir); BVP akma başlangıcı 2·Mp/s_g braketi + doyma **tam Mp** (100.0000000000) + elmas yüzeyi f=0.000e+00 (tri6+tri15; test_plate_plastic) |
| Eksenel-simetrik plate eğilmesi (§2.4) | +%0.6 | **YOK** — axisym yapısal eleman v1.x (dürüst boşluk) |
| Updated mesh / büyük deformasyon (§2.6) | eğri üstünde | **YOK** — v1.x (ROADMAP kapsam matrisi) |

## C. Doğrudan PLAXIS-sonucu kıyasları

| Problem | PLAXIS sonucu | KATAI 2D | Δ |
|---|---|---|---|
| **Tutorial Lesson 1** — kum üzerinde dairesel temel (MC + K0 + su tablosu + verilen oturma) | 588 kN | 606.9 kN | **+%3.2** |
| **Konsol palplanş, kohezif kil** (IJCRT 2024, PLAXIS Tablo 6.2; Undrained B) — iyi-konumlu vaka M_max | 19.4 kNm/m | 18.8 (tri15) | **−%3** |
| Aynı problem, drenajlı ↔ analitik LEM | LEM 14.4 | 14.6 | **+%1** |
| Derin duvarlar (D≥2 m) bandı | 26.4…66.1 | −%39…+%13 (domain ölçeğiyle PLAXIS'e yakınsar: −47→−17%) | domain dersi kayıtlı |

## D. Dinamik / sismik (PLAXIS sayı yayımlamaz → bağımsız oracle + yöntem paritesi)

| Konu | Referans | KATAI 2D sonucu |
|---|---|---|
| Katmanlı saha tepkisi (SHAKE-sınıfı lineer çekirdek) | viskoelastik transfer-matris (bağımsız; Rayleigh'in TAM sürekli karşılığı) | tek tabaka rezonans **−%0.18**; iki tabaka tepe **−%0.09** (\|T\|=17), flank +%0.06 |
| Compliant (absorbing) taban | Joyner-Chen; **PLAXIS Sci §6.3.2 + Tut §17.8.5 konvansiyonları birebir kilitli** | radyasyon oracle'ına −%0.02 (üniform), −%0.05 (iki tabaka); rijit taban 8.7× çınlama dişi |
| Free-field yan sınırlar | 1D saha tepkisi | rijit tabanda −%0.1; compliant tabanda −%4.8 (köşe örtüşmesi kayıtlı) |
| Model temel frekansı f₁ | Vs/4H + katmanlı kutup | **+%0.00 / −%0.01** (iki bağımsız yol) |
| Gerçek kayıt (El Centro 1940 NS) | yayımlı PGA 0.319 g; %5 spektrum bandı | PGA 0.31882 g birebir; tepe Sa 2.87×PGA (yayımlı 2.0–3.5 bandı) |
| Nonlineer dinamik | PLAXIS Ref §11.10.4 (plastik fazla aynı yakınsama) — yöntem paritesi | lineer-limit 1.5e-13; yarı-statik elastoplastik alan %0.12; MC histerezis yalnız akmada (Mat.Man §3.5), W_elastik/W_plastik=7e-5 |
| Tasarım spektrumları | TBDY 2018 resmî katsayılar; EC8 EN 1998-1 Tablo 3.2/3.3 | plato/köşeler kapalı-formla birebir; GUI kablolaması 0.000e+00 |

## Özet okuma

- **PLAXIS'in kendi analitiğe hatası** bu sette tipik **%0.6–6**; **KATAI aynı vakalarda %1.4–9.7**
  (aynı sınıf; iki FE'nin analitiğin iki yanına düşmesi normaldir — örn. Cox'ta PLAXIS −%2.5,
  KATAI +%3.9).
- **Doğrudan PLAXIS-sayısı kıyaslarında** KATAI −%4.0 … +%6.5 bandında; en iyi eşleşmeler
  Tutorial temeli +%3.2 ve tri15 duvar −%3.
- **Dinamikte** PLAXIS yayımlı sayı yoktur; KATAI'nin bağımsız-oracle sonuçları %0.02–0.2
  sınıfındadır ve sınır-koşulu konvansiyonları PLAXIS manuallerinden birebir kilitlidir.

## Dürüstlük notları
1. PLAXIS sayıları 2–3 anlamlı haneyle yayımlanır; % farkların ~±1 puanı bu yuvarlamanın içindedir.
2. Mesh'ler birebir aynı değildir (yayınlar mesh dosyası vermez); mesh-yakınsama yönleri her
   vakada ölçülüp not edilmiştir.
3. §3.1'de ψ=φ seçimi (slip-line kıyası) ve §2.2'de E_ref≈0 modellemesi açıkça belgelidir;
   varsayım gizlenmemiştir.
4. "YOK" satırları gerçek boşluklardır ve ROADMAP kapsam matrisinde kapanış sürümleriyle kayıtlıdır.
5. Bu tablo canlıdır: her yeni benchmark (Track 7) satır ekler; sayılar ctest regresyonlarına
   bağlı olduğundan sessizce bayatlayamaz.
