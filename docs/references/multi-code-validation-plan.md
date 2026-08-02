# Çok-Kod Doğrulama Planı — PLAXIS 2D / GEO5 / FLAC / Midas GTS NX

> ℹ️ Bu plan **`docs/internal/ROADMAP.md` Track 7 (Benchmark Programı)** tarafından yürütülür; buradaki
> felsefe (analitik altın standart, yoksa kod-bandı ortası) aynen geçerlidir. Tablodaki bazı
> "❌ YOK" girdileri tarihseldir (HSsmall, embedded beam, konsolidasyon → hepsi ✅); güncel durum
> `docs/validation/VALIDATION-SUMMARY.md`.

**Başarı kriteri (kullanıcı):** KATAI sonuçları sektörün tekelleşmiş, güvenilir 2D geoteknik
programlarıyla **uyumlu** olmalı; ideal olarak hepsinin **ortasında** bir değer → "başardık".
Tek bir kodu birebir kovalamak yerine, **bu kodların manuallerindeki örnekleri** referans alıp
KATAI'yi onların bandına oturtmak. Her özellik (insanların kontrol edeceği) ayrı doğrulanır.

**Neden bant/orta-değer?** Aynı problemde PLAXIS/GEO5/FLAC/Midas bile birbirinden farklı çıkar
(eleman tipi, mesh, integrasyon, default'lar). Analitik kapalı-form varsa o altın standart (<%1-2);
yoksa kodların bandının ortası kabul edilir. **Öncelik: matematik teori + algoritma + mimari uyumu.**

## Referans kaynaklar (manuallar, çoğu ücretsiz)
- **PLAXIS 2D** Reference / Material Models / Tutorial / Scientific Manual (Seequent/Bentley, ücretsiz PDF).
- **GEO5** Theory + Engineering Manuals (Fine, çevrimiçi; sheet-pile, MSE, bearing, settlement örnekleri).
- **FLAC / FLAC2D** (Itasca) Verification & Example problems (yayımlı doğrulama seti — analitikle kıyaslı).
- **Midas GTS NX** Analysis/Verification Manual + benchmark notları.
- **SoilVision / SVSolid, OptumG2** (limit analiz, ek çapraz-kontrol).
- Analitik/akademik: Prandtl, Reissner, Boussinesq, Lamé, Terzaghi, Rankine, Blum, Bishop/Spencer (şev).

## Özellik bazında doğrulama matrisi
İnsanların ilk bakacağı özellikler (kullanıcı listesi) + durum + referans hedefi:

| Özellik | KATAI durum | Referans benchmark (kaynak) | Hedef |
|---|---|---|---|
| **Öz-ağırlık / K0** | ✅ | jeostatik σ'v=γ'z, σ'h=K0σ'v; K0 procedure (PLAXIS Sci. Man.) | <%1 |
| **Mohr-Coulomb** | ✅ | Prandtl Nc=5.14, Reissner Nc(φ), Rankine Ka/Kp (analitik + tüm kodlar) | <%2 |
| **Hardening Soil** | ✅ büyük oranda | Berlin Sand triaxial (PLAXIS MMM Fig 15.4); oedometer Eoed | <%5 |
| **HSsmall** | ❌ YOK | küçük-şekildeğiştirme G0/γ0.7 (PLAXIS MMM §7; FLAC örnek); kazı yüzey-oturması | ekle→<%5 |
| **Yeraltısuyu (seepage)** | ✅ | confined flow-net, baraj-altı uplift (Harr analitik); freatik yüzey (Dupuit) | <%2 |
| **Dağılı yük** | ✅ | şerit yük altında Boussinesq/oturma | <%2 |
| **Nokta yük** | ⏳ kontrol | Flamant çizgi yük; tekil yük (analitik) | <%2 |
| **Plate** | ✅ (+tri15) | konsol PL³/3EI; gömülü perde sehimi/moment (PLAXIS/GEO5 sheet-pile) | <%5 |
| **Anchor** | ✅ | ankraj kuvveti (tie-back); free-earth support (GEO5 anchored wall) | <%5 |
| **Interface** | ✅ | Coulomb slip (analitik); R_inter perde sehimi etkisi (PLAXIS) | <%5 |
| **Embedded beam** | ❌ YOK | pile row (PLAXIS MMM §18.4; eksenel + yanal kapasite) | ekle |
| **Kazı (excavation)** | ✅ kuplajlı | kantilever palplanş kilde (IJCRT/PLAXIS); **tie-back kazı tutorial** (PLAXIS Lesson) | bant-orta |
| **Şev FoS** | ✅ | ACADS/Donald-Giam seti; Bishop/Spencer (GEO5, Slide2 verification) | <%5 |
| **Taşıma gücü** | ✅ | Prandtl/Reissner/Vesic Nγ (analitik + kodlar) | <%2-5 |
| **Konsolidasyon (Biot)** | ❌ YOK | Terzaghi 1D U-Tv (analitik); PLAXIS consolidation tutorial | ekle |

## Kanonik kuplajlı benchmark'lar (kazı/perde — birden çok kodda var)
1. **PLAXIS Tutorial "Excavation / tie-back wall"** — tam parametreli, çok-fazlı; perde sehimi+moment+
   ankraj kuvveti. (Çok-fazlı interface taşıma + anchor aktivasyon gerektirir → KATAI'de sıradaki altyapı.)
2. **GEO5 "Sheet Pile Wall" / "Anti-Slide Pile"** Engineering Manual — LEM + FEM, tam veri.
3. **FLAC sheet-pile / braced excavation** verification — analitik/PLAXIS kıyaslı.
4. **Schweiger / Berlin excavation benchmark** (akademik, çok-kod kıyas — HSsmall ile).

## Doğrulama disiplini (her benchmark için)
1. Problemi TAM ve adım adım doğru kur (geometri, BC, malzeme, mesh, fazlar) — kurulum hatası = en sık hata
   kaynağı (bkz domain-boyutu dersi: sınırlar kırılma mekanizmasını kesmemeli).
2. Mesh + domain yakınsama kontrolü (sonuç sayısaldan değil fizikten gelmeli).
3. Mümkünse aynı eleman (tri15) + aynı malzeme modeli + aynı mesh felsefesi (graded) = sektör uyumu.
4. Sonucu kaynak bandıyla kıyasla; analitik varsa <%1-2, yoksa kod-bandının ortası.
5. `docs/validation/<benchmark>.md`'ye tam parametre + sonuç + değerlendirme yaz; sayısal-cevaplılar
   regresyon paketine.

İlgili: [[validation-benchmarks]], [[material-model-architecture]], [[literature-review]],
`docs/validation/wall-benchmark-plaxis.md`.
