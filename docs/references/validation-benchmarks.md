# Doğrulama Benchmark'ları — Referanslar

Hedef: kritik sonuçlarda PLAXIS'ten **< %5** sapma. Sıra: tek eleman → analitik
süreklilik → plastisite limit yükü → kuplajlı → devlerle birebir.

## Seviye 0 — Tek eleman (bünye modeli)
- Bilinen gerilme yolu (triaxial, oedometer) uygula; return-mapping çıktısını
  elle/Excel ile karşılaştır. MC ve HS hatalarını mesh'e girmeden yakalar.
- HS için: drained triaxial (E50 doğrulama), oedometer (Eoed), un/reloading (Eur).

## Seviye 1 — Lineer elastik analitik (kapalı form)
- **Boussinesq:** yarım-uzayda tekil/şerit yük → düşey gerilme, oturma.
- **Flamant:** çizgi yük.
- **Lamé (kalın silindir):** **axisymmetry** doğrulaması (radyal/teğet gerilme).
- **Kirsch (delikli plaka):** gerilme yoğunlaşması (kσ=3) → gerilme doğruluğu.
- **Konsol kiriş uç sehimi:** plate/beam yapısal elemanı.
Kaynak: herhangi elastisite kitabı (Timoshenko & Goodier), P&Z-1.

## Seviye 2 — Plastisite limit yükleri (MC)
- **Prandtl (1921):** ağırlıksız zeminde şerit temel → **Nc = 2 + π = 5.14**
  (pürüzlü taban). Plastisitenin en keskin testi.
- **Terzaghi/Prandtl taşıma gücü:** Nc, Nq, Nγ faktörleri. (Nq = e^(π·tanφ)·tan²(45+φ/2)).
- **Rankine/Coulomb:** aktif/pasif toprak basıncı katsayıları → perde.
- **Sonsuz şev:** analitik FoS = (c + γz cos²β tanφ)/(γz sinβ cosβ).
Kaynaklar:
- Bearing capacity by FE (Sloan/Newcastle): https://www.newcastle.edu.au/__data/assets/pdf_file/0006/22596/72_Numerical-limit-analysis-solutions-for-the-bearing-capacity-factor-N-gamma.pdf
- Nc=5.14 FE doğrulama (NS-FEM): https://link.springer.com/chapter/10.1007/978-981-15-2184-3_147

## Seviye 3 — Klasik geoteknik benchmark
- **Terzaghi 1D konsolidasyon:** U–Tv ilişkisi (analitik, Fourier serisi).
  Pratik: Tv = (π/4)(U)² (U≤%60); Tv = −0.933·log(1−U) − 0.085 (U>%60).
  Tv = cv·t / d² (d = drenaj yolu). Kaynak:
  https://scirp.org/html/8-1880191_49800.htm
- **Baraj/palplanş altı sızma:** kararlı akış, flow-net, sızma debisi.
- **ACADS şev benchmark seti (1989, Donald & Giam):** 5 temel + 5 varyant problem;
  limit denge (Bishop/Spencer) referans FoS'larıyla karşılaştırma. Rocscience
  Slide2 Verification Manual bunları çözülmüş cevaplarıyla içerir →
  https://static.rocscience.cloud/assets/verification-and-theory/Slide2/Slide_SlopeStabilityVerification.pdf
- **Konsol/ankrajlı palplanş:** Blum analitik + PLAXIS dersi.
- **KANTİLEVER PALPLANŞ — KİLDE (IJCRT2024, Paul/Halder/Mukherjee) — KİLİTLİ PLAXIS BENCHMARK:**
  Saf kohezyonlu kil (φ=0, Undrained B), kantilever palplanş, tek-aşama kazı, SUSUZ.
  *Zemin:* γ=17 kN/m³, E'=150 MPa, ν'=0.4, Cu=25/30/35 kPa, φ=0, ψ=0.
  *Palplanş (PU-12-240, Arcelor Mittal, elastik):* EA=2.94×10⁶ kN/m, EI=45360 kN·m²/m, w=1.101, ν=0.28.
  *Interface:* R_inter=0.67. z0=2Cu/γ (tension-crack derinliği). Gömme D = LEM'den (×1.4 FoS).
  *PLAXIS 2D bending moment (Tablo 6.2) — STABİL vakalar (H[m]/Cu[kPa]→D[m], BM[kNm]):*
  4.0/25→1.33, 7.9 · 4.5/25→3.29, 26.4 · 4.5/30→0.91, 3.5 · 5.0/30→2.17, 19.4 · 5.5/30→4.52, 44.7 ·
  5.5/35→1.57, 12.5 · 6.0/35→3.11, 42.1 · 6.5/35→5.80, 66.1. (LEM Tablo 5.1: aynı H/Cu→ örn 5.0/25 BM=70.7.)
  "--" vakaları (3.0/25, 5.0/25, 6.0/30, 7.0/35) = kritik yükseklik aşıldı/çöktü (PLAXIS yakınsamadı).
  NOT: K0 belirtilmemiş (φ=0 → Jaky K0=1 veya ν/(1−ν)=0.667); net basınç pa=σv−2c/pp=σv+2c φ=0'da K0'a
  görece duyarsız. Kaynak: IJCRT Vol12 Issue7 (2024) "Analysis Of Cantilever Sheet Pile Embeded In Cohesive
  Soil", https://www.ijcrt.org/papers/IJCRT21X0271.pdf. (LEM hand-calc'ları özensiz/tutarsız γ; PLAXIS kurulumu
  Ch6 eksiksiz → FE-vs-FE kıyas hedefi PLAXIS BM Tablo 6.2.)

## Seviye 4 — Ticari/akademik kodlarla birebir
- **PLAXIS 2D tutorial dersleri** (temel oturması, kazı, dolgu) → birebir tekrar.
  PLAXIS Reference/Tutorial Manual (Seequent, ücretsiz).
- **OpenSees** (BSD): zemin kolonu, element-seviyesi bünye doğrulaması. Sonuç
  karşılaştırması lisans açısından serbest.
- **GEO5 / Midas GTS NX:** nihai çapraz-kontrol (lisans varsa).
- **SGM** çalışan kodlarının çıktıları (elastisite, plastisite, konsolidasyon):
  https://inside.mines.edu/~vgriffit/PFEM5

## Regresyon paketi
Yukarıdakilerin sayısal-cevaplı olanları otomatik test paketine girer; her build'de
koşar, bir sonuç %5 toleransı aşarsa build **kırmızı**. (bkz STACK.md — CI/build.)
