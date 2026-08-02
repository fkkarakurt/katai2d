# Gerçek Deprem Kaydı — Akselerogram Girişi + El Centro 1940 NS (Benchmark Dalga-1)

**Track 7 / ROADMAP v0.6 kapısı.** Test: `test_real_record` (ctest). Tarih: 2026-07-20.
Risk sınıfı: **Vital** (gerçek tasarım koşuları artık gerçek kayıtla sürülür).

## Özellik: akselerogram girişi (SeismicWave::Record)

Dinamik faz artık kullanıcı akselerogramı alır: `Phase.accel_record` [m/s²] + `record_dt`;
`seismic_amp` ölçek çarpanı (1 = kayıt aynen); kayıt bitince girdi SIFIR (serbest titreşim —
dürüst, sarılmış/tutulmuş sinyal değil). GUI: dalga combo'sunda "Accelerogram (record)" + Import
(iki-sütun t/a veya satır-başı-ivme; birim oto-algılama: |a|max < 2 → g kabul edilip 9.81 ile
çarpılır — gerçek bir kaydın m/s² tepesi 2'nin altına, g tepesi üstüne çıkmaz; içe aktarım sonrası
PGA iki birimde de gösterilir = yanlış tahmin bir bakışta görünür). Kayıt PROJE DOSYASINDA saklanır
(koşu, başıboş yan dosyasız yeniden-üretilebilir). Adım uyarısı kayda özel: çözücü dt'si kayıt
dt'sinden kabaysa "kayıt ALT-ÖRNEKLENİYOR" uyarısı (gerçek kriter, uydurma eşik değil).

## V&V

**(a) Kablolama kimliği:** örnekleri tam A·sin(2πf·k·dt) olan bir Record koşusu, Harmonic koşuyu
**1.9e-12 bağıl** farkla üretir (yuvarlama içi — (step·dt)/dt'nin ~1 ulp kayması + sin argüman
birleşim sırası; 1e-9 üzeri gerçek kablolama hatası olurdu, bant dürüstçe budur).

**(b) Kayıt kimliği:** depodaki `tests/data/elcentro-1940-ns.dat` (köken: yanındaki .md) 1560 örnek,
dt=0.02 s, t=0→31.18 s; **PGA = 0.31882 g @ t=2.02 s** — yayımlanmış klasik 0.319 g @ ~2 s ile
birebir (sayısallaştırmalar arası ±%2-3 doğaldır; bu dosyanın rolü "kimliği doğrulanmış gerçek
kayıt" olmaktır).

**(c) %5 sönümlü tepki spektrumu (doğrulanmış response_spectrum motoruyla):** tepe Sa = **0.916 g
@ T=0.19 s = 2.87×PGA** — yayımlanmış %5 El Centro büyütme bandının (2.0–3.5×PGA) içinde; tepe
kısa-orta periyot bölgesinde (yayımlanmış grafiklerde ~0.15–0.7 s arası karşılaştırılabilir yerel
maksimumlar; hangisinin global çıktığı sayısallaştırmaya göre değişir — dürüst bant); uzun-periyot
ordinatı küçük (Sa(3 s)=0.118 g < 0.15 g). Nokta-nokta yayımlanmış tablo kıyası AÇIK kalem olarak
kayıtlıdır (motorun kendisi kapalı-formlarla ayrıca pinli: seismic-verification.md).

**(d) Ürün koşusu (iki-tabaka benchmark profili, gerçek kayıt, 1560 adım):**
rijit taban tepe yüzey ivmesi **1.59 g**, compliant (yutucu) taban **0.52 g** — radyasyonun gerçek
kayıttaki etkisi; yutucu taban tam-yansıtıcıyı aşmaz (fiziksel yön), iki koşu da spektrum üretir.

## Reprodüksiyon
`cmake --build build/msvc-rwdi --target test_real_record` → çalıştır (veri: tests/data, CMake
`KATAI_TEST_DATA_DIR` ile bağlanır). Kaynaklar: vibrationdata.com/elcentro.dat (veri) +
klasik literatür PGA/spektrum karakteristikleri (Kramer 1996; Chopra).
