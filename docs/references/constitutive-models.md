# Bünye Modelleri (Constitutive Models) — Referanslar

Öncelik sırası: **Linear-Elastic → Mohr-Coulomb → Hardening Soil**.

## 1. Linear-Elastic (LE)
- Klasik Hooke; iki parametre (E, ν) veya (G, K). Plane strain & axisymmetry için
  D matrisleri standart — bkz ZT, SGM Ch 5.
- Doğrulama: analitik elastisite (Boussinesq, Lamé, Kirsch) — validation dosyasında.

## 2. Mohr-Coulomb (MC) — elasto-plastik, mükemmel plastik
**Formülasyon kaynağı:** PLAXIS-MM (MC bölümü), P&Z-1, dSN.

**Akma yüzeyi:** altı düzlemli piramit (principal stress uzayında köşeli).
Parametreler: c (kohezyon), φ (içsel sürtünme), ψ (dilatasyon), E, ν.
**Non-asosiye akış** (ψ ≠ φ) standarttır.

**Asıl zorluk — köşe ve apeks tekilliği:** akma yüzeyinin kenarlarında ve tepesinde
gradyan süreksiz; stres entegrasyonu burada patlar. İki ana strateji:
- **Clausen ve ark. — principal stress uzayında return mapping.** Regular / kenar
  (sağ-sol) / apeks dönüşlerini ayrı ayrı, indikatörle ele alır. En temiz ve hızlı
  yaklaşım. → "An Efficient Return Algorithm for Non-Associated Mohr-Coulomb
  Plasticity" (Clausen, Damkilde, Andersen).
- **Abbo & Sloan — C2-sürekli (hiperbolik/yuvarlatılmış) MC yüzeyi.** Apeksi
  yuvarlayıp tekilliği tamamen ortadan kaldırır; Newton-Raphson çok daha kararlı.
  → Newcastle (Sloan grubu) PDF'leri.
- Alternatif: subdifferential tabanlı implicit return mapping (arXiv 1508.07435).

**Implementasyon kararı (2026-06-03, rafine):** **principal-stress uzayında
kapalı-form return mapping** (Clausen et al. 2007 / Sysala & Čermák 2016). Abbo-Sloan
yuvarlatması elendi (yaklaşım parametresi getiriyor; "PLAXIS-üstü doğruluk" için
yaklaşımsız form yeğ). Mükemmel-plastik (H=0) ilk sürümde Δλ afin → iterasyonsuz.
Dört bölge (düz yüzey/sol kenar/sağ kenar/apeks) ve tutarlı teğet **doğrulanmış
denklemlerle** `mohr-coulomb-formulation.md` içinde sabitlendi. Doğrulama: material-point
(triaxial, f≈0'a tam dönüş) → Prandtl Nc=5.14.

Kaynaklar:
- Clausen et al., return algorithm: https://www.researchgate.net/publication/269155657
- Abbo & Sloan, C2 MC yüzeyi: https://www.newcastle.edu.au/__data/assets/pdf_file/0013/22252/A-C2-continuous-approximation-to-the-Mohr-Coulomb-yield-surface.pdf
- Return mapping (Sloan grubu özet): https://www.newcastle.edu.au/__data/assets/pdf_file/0016/22309/Return-Mapping-Algorithms-and-Stress-Predictors-for-Failure-Analysis-in-Geomechanics.pdf
- Subdifferential implicit return mapping: https://arxiv.org/pdf/1508.07435

## 3. Hardening Soil (HS) — projenin en zor parçası
**Birincil kaynak:** Schanz, Vermeer & Bonnier (1999), "The Hardening Soil Model:
Formulation and Verification", *Beyond 2000 in Computational Geotechnics*, Balkema,
s. 281–296. + PLAXIS-MM (HS bölümü).

**Çekirdek fikir:**
- Gerilmeye bağlı rijitlik (stress-dependent stiffness): E50, Eoed, Eur ayrı; güç
  yasası `m` ile (`E = Eref·((c·cotφ + σ)/(c·cotφ + pref))^m`).
- **Çift sertleşme (double hardening, Vermeer 1978):** deviatorik (kayma) akma
  yüzeyi + hacimsel (cap) akma yüzeyi. İzotropik sertleşme, plastik kayma + hacim
  şekildeğiştirmesine bağlı.
- Deviatorik kısım Duncan-Chang hiperbolik ilişkisini genelleştirir; cap için
  asosiye akış + dilatasyon.

**Uyarı (kararımız D8):** PLAXIS'in tam implementasyonu kısmen tescilli. Biz
yayınlanmış denklemlere dayanacağız → sonuçlar PLAXIS'ten *biraz* sapabilir. HS en
çok doğrulama gerektiren modeldir; tek-eleman testleriyle (triaxial, oedometer
simülasyonu) başlayacağız.

**Sonraki adım (ileride):** HS-small (HSsmall) — Benz (2007) doktora tezi,
"Small-Strain Stiffness of Soils and its Numerical Consequences", Stuttgart.

Kaynaklar:
- Schanz-Vermeer-Bonnier (geotechpedia): https://geotechpedia.com/Publication/Show/975/THE-HARDENING-SOIL-MODEL--FORMULATION-AND-VERIFICATION
- Semantic Scholar (PDF): https://www.semanticscholar.org/paper/ced1b3f33e91eeb83951913262611e5a3e85bf6b
- Geoengineer eğitim notu (HS özet): https://www.geoengineer.org/education/numerical-constitutive-modeling/numerical-modelling-the-hardening-soil-model

## 4. Undrained davranış (PLAXIS yöntemleri A/B/C)
**Kaynak:** PLAXIS-MM (Undrained behaviour bölümü).
- **Undrained (A):** efektif gerilme analizi, **efektif** mukavemet (c', φ'). Su
  basıncı modelden çıkar. Önerilen yöntem; konsolidasyonla mukavemet artışını yakalar.
- **Undrained (B):** efektif gerilme analizi, **undrained** mukavemet (su, φ=0).
- **Undrained (C):** total gerilme analizi, undrained parametreler (Eu, νu≈0.495, su).
- Mekanizma: efektif G, φ' → undrained Eu, νu'ya dönüştürülür; su için yüksek
  bulk modülü (Kw/n) eklenir. **Hacimsel kilitlenmeyi (volumetric locking)** önlemek
  için yüksek mertebeli eleman (15-düğümlü) önemli.

Kaynak (resmi manuel): PLAXIS-MM (yukarıdaki Seequent PDF).
