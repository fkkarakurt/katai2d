# Elemanlar ve Mesh Üretimi — Referanslar

## 1. Eleman kütüphanesi
**Karar (D10/Q3):** kullanıcı seçer — **6-düğümlü** (linear strain) ve **15-düğümlü**
(4. mertebe) üçgen. PLAXIS varsayılanı 15-düğümlü.

- **15-düğümlü üçgen:** 12 Gauss noktası; PLAXIS'in yüksek doğruluk + undrained/
  incompressible problemlerde **hacimsel kilitlenmeye direnç** için tercih ettiği
  eleman. Bizim <%5 hedefimizin anahtarı.
- **6-düğümlü üçgen:** 3 Gauss noktası; daha hızlı, basit problemler için yeterli.
  Kilitlenmeye daha yatkın.
- Yapısal elemanlar mesh ile **uyumlu** olmalı: 15-düğümlü üçgenin kenarında
  5-düğümlü hat elemanı (plate/geogrid), 6-düğümlünün kenarında 3-düğümlü.
- Şekil fonksiyonları, Gauss kuralları, plane-strain/axisymmetry B & D matrisleri:
  ZT, SGM Ch 4–6, P&Z-1 Ch 2.

**Genel FEM teorisi:** Zienkiewicz & Taylor; Bathe, *Finite Element Procedures*.
IFEM ders notları (üçgen elemanlar, açık erişim):
https://quickfem.com/wp-content/uploads/IFEM.Ch15.pdf

## 2. Mesh üretimi (sıfırdan yazıyoruz — D5 istisnası)
Permissive + kaliteli hazır kütüphane yok (Triangle: non-commercial; Gmsh/CGAL:
GPL). Bu yüzden **kendi mesher'ımızı** yazacağız. Algoritma referansı:

- **Shewchuk — "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay
  Triangulator".** *Kodu kullanmıyoruz* (lisans), ama bu makale + tez kısıtlı
  Delaunay + Ruppert iyileştirmesinin en iyi mühendislik anlatımı.
  PDF: https://people.eecs.berkeley.edu/~jrs/papers/triangle.pdf
- **Ruppert's Delaunay Refinement:** PSLG (planar straight-line graph) al, sadece
  kaliteli üçgen üret (circumradius/en-kısa-kenar oranı eşiği). Quake sayfası:
  https://www.cs.cmu.edu/~quake/tripaper/triangle3.html
- **Chew's 2. algoritması:** kısıtlı Delaunay üzerinde min açı ~28.6°'ye kadar
  garanti. Alternatif iyileştirme stratejisi.
- Shewchuk, "Delaunay Refinement Algorithms for Triangular Mesh Generation":
  https://people.eecs.berkeley.edu/~jrs/papers/2dj.pdf
- Wikipedia özet (hızlı bakış): https://en.wikipedia.org/wiki/Delaunay_refinement

**Yapı taşları (kendi implementasyonumuz):**
1. Kısıtlı Delaunay triangülasyon (zemin tabaka sınırları, yapı kenarları = kısıt).
2. Ruppert/Chew kalite iyileştirmesi (min açı, max alan; bölgesel yoğunlaştırma).
3. 6→15 düğüm zenginleştirme (ara düğüm ekleme, kenar uyumu).
4. Mesh kalite metrikleri + iyileştirme (Laplacian/optimization smoothing).

**Quad mesh (ikincil, sonra):** düzgün bölgeler için mapped/structured quad;
keyfi geometride paving. MVP'yi bloklamaz.

**Geometri çekirdeği (permissive kütüphaneler):**
- Poligon boolean (tabaka kesişimleri, kazı çıkarma): **Clipper2** (Boost lisans).
- Geometrik yardımcılar: **Boost.Geometry** (Boost lisans).
