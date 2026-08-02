# Analiz Tipleri ve Yapısal Elemanlar — Referanslar

## 1. Nonlineer çözüm (tüm plastisite/aşamalı analizlerin altyapısı)
- **Newton-Raphson / modified N-R**, yük/deplasman kontrolü; yakınsama
  kriterleri. Kararsız (softening, limit yük) durumlarda **arc-length** (Crisfield).
- Kaynak: ZT, P&Z-1 Ch 9, SGM Ch 6, Simo-Hughes. dSN (pratik return mapping +
  consistent tangent operatörü → kuadratik yakınsama).

## 2. Strength Reduction (φ-c reduction) — güvenlik/şev
**Kanonik kaynak:** Griffiths & Lane (1999), "Slope stability analysis by finite
elements", *Géotechnique* 49(3):387–403. (SGM yazarlarından Griffiths.)
- c ve tanφ aynı anda SRF ile bölünür; **yakınsama bozulana** kadar artırılır →
  SRF = güvenlik sayısı. Kayma yüzeyini önceden varsaymaya gerek yok (limit denge
  yöntemine göre avantaj).
- φ<30° / küçük sürtünme açılarında modifiye yaklaşımlar (MSRFEM) var.
Kaynaklar:
- ScienceDirect (verimli SRM): https://www.sciencedirect.com/science/article/pii/S0266352X24005329
- Springer (yükseklik/açı etkisi): https://link.springer.com/article/10.1186/s43088-021-00115-w

## 3. Konsolidasyon (kuplajlı, zamana bağlı)
- **Biot (1941)** kuplajlı teori: deplasman–boşluk suyu basıncı (u–p) formülasyonu.
  PDE seti: elastisite + Darcy süreklilik.
- **Karışık eleman / inf-sup (LBB) uyarısı:** eşit-mertebeli u ve p alanları
  kararsızdır (basınç salınımı). Çözüm: deplasman için yüksek mertebe (15-düğüm),
  basınç için daha düşük mertebe; veya stabilize formülasyon. → 15-düğümlü
  elemanda 6 düğümde pore pressure standart yaklaşımdır (PLAXIS de böyle yapar).
- Kaynak: P&Z-1 Ch 10 (kuplajlı), SGM Ch 9, Biot 1941.
- Analitik karşılaştırma: Terzaghi 1D konsolidasyon (validation dosyası).
Kaynak (mixed FE özet): https://www.researchgate.net/publication/222660034

## 4. Yeraltı suyu akışı (sızma)
- **Kararlı akış (steady-state):** Darcy + Laplace; kapalı (confined) akış basit.
- **Serbest yüzey (unconfined) / freatik yüzey:** ana zorluk — freatik yüzey
  konumu baştan bilinmez, **iteratif** bulunur (ıslak/kuru bölge sınırı). Sabit-grid
  yöntemleri (permeabiliteyi basınca göre güncelle) keyfi geometride pratiktir.
- Kaynak: P&Z (akış), SGM (akış bölümü).
- Serbest yüzey FE (pratik yöntem): https://www.sciencedirect.com/science/article/abs/pii/S0266352X02000034

## 5. Interface elemanları (soil-structure)
**Kaynak:** Goodman, Taylor & Brekke (1968) — zero-thickness joint element;
Day & Potts (1994), "Zero thickness interface elements—numerical stability and
application", *IJNAMG* 18:689–708.
- Düğüm çiftleri **çakışık** (zero-thickness); kayma + ayrılma (sliding/separation/
  debonding). Normal & kayma rijitliği (kn, ks) → **suni rijitlik ölçekleme**
  sorununa dikkat (çok büyük kn sayısal sorun çıkarır).
- PLAXIS: bitişik zemine uygulanan **Rinter** (mukavemet/rijitlik azaltma faktörü).
- Day & Potts: çakışık düğümlü interface'in sayısal kararlılık tuzakları.
Kaynaklar:
- Day & Potts (1994): https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.1610181003
- Zero-thickness interface (özet): https://www.researchgate.net/publication/263298303

## 6. Yapısal elemanlar
- **Plate (perde/kiriş):** Mindlin kiriş (kayma deformasyonlu), mesh ile uyumlu hat
  elemanı. EA, EI, w. Kaynak: P&Z-2, PLAXIS Reference Manual.
- **Geogrid:** sadece çekme (tension-only) taşıyan eksenel eleman; EA.
- **Anchor:** node-to-node (elastoplastik yay) ve fixed-end ankraj; EA, Fmax.
- **Embedded pile row (2D'de kazık):** Sluis & Brinkgreve ve ark. (2014); Sluis
  MSc tezi (Delft, 2012) "Validation of Embedded Pile Row in Plaxis 2D". Kiriş +
  özel **interaksiyon arayüzü** (mesh'e gömülü); 2D'de düzlem-dışı kazık aralığı
  (Lspacing) ile çalışır. Turello ve ark. (2016) interaksiyon yüzeyi iyileştirmesi.
Kaynaklar:
- Embedded pile row (Plaxis 2D, Sluis): https://bouassidageotechnics.wordpress.com/wp-content/uploads/2016/05/embedded-pile-row-plaxis.pdf
- Validation (ResearchGate): https://www.researchgate.net/publication/283295845
- Improved embedded beam (Smulders, Hosseini, Brinkgreve 2019): https://www.issmge.org/uploads/publications/51/75/0193-ecsmge-2019_Smulders-Hosseini-Brinkgreve.pdf
