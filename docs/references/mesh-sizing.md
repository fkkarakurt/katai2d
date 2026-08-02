# Yerel Mesh Yoğunluğu (Sizing Field / Coarseness) — Formülasyon

PLAXIS 2D'nin yerel mesh kontrolü semantiği: her geometrik nesnenin (bölge/poligon, yapısal hat,
yük) bir **Coarseness factor** f'si vardır (varsayılan 1; küçük = ince). KATAI bunu Ruppert
iyileştiricisine **konuma bağlı alan sınırı** (sizing field) olarak verir.

Kaynaklar: Shewchuk (2002) *Delaunay Refinement Algorithms for Triangular Mesh Generation*
(alan/boyut kısıtlı Ruppert, sonlanma); Persson (2006) *Mesh size functions* (Lipschitz-gradasyonlu
boyut alanı); PLAXIS 2D Reference Manual — Mesh (Coarseness factor, refine/coarsen, yapı/yük
çevresinde otomatik sıklaştırma). PLAXIS'in içsel sabitleri kamuya kapalı; buradaki sabitler
KATAI'nin kendi (dokümante) seçimleri, semantik PLAXIS ile eş.

## 1. Boyut alanı
Hedef kenar boyu `h0 = sqrt(2·max_area)` (GUI eleman boyu). Nesne katsayıları `f ∈ [1/16, 4]`
(clamp; pozitif alt sınır ⇒ Ruppert sonlanması garanti):

```
h_region(x) = h0 · f_poly(x)            (içinde bulunduğu son poligon kazanır — malzeme gibi)
h_src,i     = h0 · f_i · (auto ? 0.5 : 1)   (yapısal hat / yük çizgisi / tekil yük kaynağı)
h(x)        = min( h_region(x), min_i [ h_src,i + g · dist(x, kaynak_i) ] )      g = 0.5
alan(x)     = h(x)² / 2                  (üçgenin CENTROID'inde değerlendirilir)
```

- **Gradasyon g = 0.5:** kaynaktan uzaklaştıkça boyut en çok 0.5 m/m büyür (Lipschitz);
  komşu elemanlar arasında ani boyut sıçraması olmaz, geçiş bandı pürüzsüz.
- **Auto-refine (EMR-vari, GUI varsayılanı AÇIK):** plate/geogrid/embedded-beam hatları ve
  yükler kendi f'lerine ek ×0.5 alır (PLAXIS'in yapı çevresi otomatik sıklaştırması gibi).
  API varsayılanı KAPALI (mevcut test mesh'leri bit-birebir korunur); GUI açık geçer.
- **Refine/Coarsen (sağ-tık):** f ×= 0.5 / f ×= 2 (PLAXIS komut semantiği), anında re-mesh.

## 2. Mesher entegrasyonu
`quality_mesh(..., SizeField, ...)`: `find_bad_triangle` alan sınırını sabit `max_area` yerine
centroid'de `max_area_at(gx, gy)` ile değerlendirir; encroachment/min-açı/sınır mantığı birebir
aynı. Alan boş (factor'sız) modelde hiç kurulmaz → eski yol **bit-birebir** (test A).

## 3. Doğrulama (`test_mesh_refine`)
| Kontrol | Sonuç |
|---|---|
| (A) factor'sız = eski sabit-alan mesh | eleman/düğüm sayısı birebir |
| (B) bölge f=0.25 | bölge içi max alan = cap (0.1250/0.1250); üst bölge iri kalır |
| (C) plate auto-refine ×0.5 | hat dibinde cap'e uygun (0.55 ≤ 0.845 gradasyonlu sınır), uzakta 1.89'a büyür |
| (D) Ruppert min-açı | 22.1° ≥ 20° (değişken alan kaliteyi bozmaz) |
| (E) domain bütünlüğü | hiçbir centroid dışarıda değil |

GUI: Mesh panelinde min-açı + ortalama kalite (q = 4√3·A/Σl², eşkenar=1) istatistikleri.
