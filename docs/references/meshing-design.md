# Mesher Tasarımı — Doğrulanmış Mimari ve Karar Kaydı (P1.4)

KATAI 2D mesher'ı sıfırdan yazılıyor (D5 istisnası: permissive + kaliteli hazır
kütüphane yok — Triangle non-commercial, Gmsh/CGAL GPL). Bu, PLAXIS'i **geçebileceğimiz**
sayılı kategorilerden biri: PLAXIS 2D eski bir mesh kütüphanesi kullanıyor. Üç kaldıraç:

1. **Sağlamlık (robustness):** *exact/adaptive geometric predicates* (Shewchuk). Eski
   kütüphaneler kayan-nokta epsilon hack'leriyle dejenere girdide bozuk/ters eleman
   üretir veya sonsuz döngüye girer. Exact predicates ile **asla**.
2. **Kalite:** Ruppert/Chew Delaunay refinement — **ispatlı** minimum açı sınırı
   (sliver yok) → daha iyi matris koşullanması → daha doğru FEM + daha hızlı çözüm.
3. **Hız:** filtered predicates (hızlı yol), BRIO/Hilbert sıralı incremental ekleme
   (beklenen ~lineer), üçgen-tabanlı komşuluk + veri-odaklı SoA (cache-dostu).

## Kilitlenen kararlar
- **Geometric predicates:** Shewchuk *adaptive* `orient2d` + `incircle` — statik hata
  filtresi (hızlı yol) + exact expansion aritmetiği (geri-düşüş, `std::fma` tabanlı
  TwoProduct/TwoSum). Doğruluk = işaret kesinliği; dejenere (collinear/cocircular)
  girdide tam 0.
- **Triangülasyon kurulumu:** **incremental Bowyer-Watson** + uzamsal sıralama
  (BRIO/Hilbert). Gerekçe: Steiner noktası eklemeyi (refinement) doğal destekler →
  ilk triangülasyon ve iyileştirme **aynı** çekirdek. (Divide-and-conquer elendi:
  statik nokta kümesi için hızlı ama artımlı eklemeyi taşımıyor.)
- **Kısıtlı Delaunay (CDT):** segment ekleme + kenar çevirme (Sloan 1993 / Shewchuk).
  Zemin tabaka sınırları, yapı kenarları, kazı hatları = kısıt segmentleri.
- **Kalite iyileştirme:** **Ruppert** (encroachment → segment bölme; skinny üçgen →
  çevrel merkez ekleme), *local feature size*'a göre kademeli (graded). Hedef min açı
  ~ ≤ ~30° (Ruppert ispatlı ~20.7°, pratikte daha yüksek; Chew 2. ~28.6° alternatif).
- **Veri yapısı:** üçgen-tabanlı, üçgen başına 3 köşe + 3 komşu üçgen (O(1) gezinme),
  SoA-dostu. Nokta-konumlandırma: komşu-yürüme (walk).
- **Boyutlandırma (sizing):** geoteknik ihtiyaç — yük/kazı köşelerinde yoğunlaştırma;
  kullanıcı boyut fonksiyonu + local feature size.
- **Geometri yardımcıları (permissive):** poligon boolean (tabaka kesişimi, kazı
  çıkarma) için **Clipper2** (Boost lisans) ileride; şimdilik kendi PSLG'miz.

## Alt-adımlar (her biri ayrı doğrulanır)
- **P1.4a — Sağlam predicates:** `orient2d` (+ `incircle`). Doğrulama: işaret kesinliği,
  collinear/cocircular → tam 0, dejenereye-yakın girdide naive double yanılırken doğru
  işaret, antisimetri.
- **P1.4b — Incremental Delaunay çekirdeği:** nokta kümesi → Delaunay (Bowyer-Watson +
  uzamsal sıralama). Doğrulama: Delaunay özelliği (boş çevrel-çember), Euler ilişkisi,
  konveks örtü, rastgele büyük kümede sağlamlık.
- **P1.4c — Kısıtlı Delaunay:** segment ekleme/kurtarma. Doğrulama: kısıtların korunması,
  CDT özelliği.
- **P1.4d — Ruppert iyileştirme + boyutlandırma:** kalite (min açı, max alan), kademeli.
  Doğrulama: açı histogramı sınırda, sonlanma, boyut-optimumluk.
- **P1.4e — tri6/tri15 düğüm üretimi:** köşe + orta-kenar (ve tri15) düğümleri, sınır
  kümeleri, malzeme bölgeleri → mevcut `Mesh` yapısına.

## Kaynaklar (doğrulama temeli)
- Shewchuk, *"Triangle: Engineering a 2D Quality Mesh Generator and Delaunay
  Triangulator"* — https://people.eecs.berkeley.edu/~jrs/papers/triangle.pdf
- Shewchuk, *"Delaunay Refinement Algorithms for Triangular Mesh Generation"* —
  https://people.eecs.berkeley.edu/~jrs/papers/2dj.pdf
- Shewchuk, *"Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
  Predicates"* (orient2d/incircle adaptive) — predicates temeli.
- Ruppert (1995), *"A Delaunay Refinement Algorithm for Quality 2D Mesh Generation"*,
  J. Algorithms. Chew, 2. algoritma (min açı ~28.6°).
- Sloan (1993), kısıtlı Delaunay triangülasyon (segment kurtarma).
