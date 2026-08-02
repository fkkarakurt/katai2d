# Tasarım kodu uyumu — Eurocode 7 (EN 1997-1) + TBDY 2018 (referans kilidi)

> v0.3 iş kolu **B** kaynağı. Katsayı değerleri **koda gömülü sabit değil, veri (config)** olacak;
> bu belge yetkili değerleri kilitler. **Doğrulanmamış hiçbir değer koda authoritative girmez.**
> Durum işaretleri: ✅ doğrulandı (çoklu kaynak) · ⏳ birincil-kaynaktan-pinlenecek.

## 0. Temel mühendislik kararı — FEM'e eşleme

İki katsayı-uygulama paradigması vardır ve FEM'e **farklı** oturur:

- **Malzeme-katsayılı (strength factoring)** — dayanım parametrelerine katsayı: `tanφ'_d = tanφ'/γ_φ'`,
  `c'_d = c'/γ_c'`, `c_u,d = c_u/γ_cu`. FEM'de **doğrudan** yapılır: azaltılmış dayanımla yeniden çöz
  (mevcut `strength_reduction.hpp` makinesi — sabit katsayılı SRM). → **EC7 DA1-C2, DA3**.
- **Direnç-katsayılı (resistance factoring)** — hesaplanan **dirence** katsayı: `R_d = R_k/γ_R`.
  Genel deformasyon FEM'i tek dirençli bir mekanizma vermez; bu **rapor/kontrol katmanında** yapılır
  (FEM çıktısı → mekanizma direnci → `/γ_R` → `E_d ≤ R_d`). → **EC7 DA2 ve TBDY 2018**.

**Sonuç:** B1a (çekirdek, FEM-native) = EC7 DA1-C2 & DA3 malzeme-katsayısı; B1b (rapor/kontrol) =
EC7 DA2 + TBDY 2018 `E_d ≤ R_d`. Bu belge her ikisinin de katsayılarını kilitler.

## 1. EN 1997-1 Annex A — önerilen (recommended) kısmi katsayılar ✅

Kaynaklar: EN 1997-1:2004+A1:2013 Annex A; UK NA to EN 1997-1; Designers' Guide to EN 1997-1;
çoklu ikincil kaynak (civils.ai, deepexcavation, MIDAS EC7 notu) — değerler çapraz-doğrulandı.
**Not:** bunlar EN önerilen değerleridir; Ulusal Ek farklılaştırabilir (veri modeli buna hazır olacak).

### 1.1 Etkiler (Actions) — Tablo A.3 ✅
| Katsayı | A1 | A2 |
|---|---|---|
| γ_G kalıcı, elverişsiz (unfav) | **1.35** | **1.00** |
| γ_G kalıcı, elverişli (fav) | 1.00 | 1.00 |
| γ_Q hareketli, elverişsiz | **1.50** | **1.30** |
| γ_Q hareketli, elverişli | 0 | 0 |

### 1.2 Zemin/malzeme parametreleri — Tablo A.4 ✅
| Katsayı | M1 | M2 |
|---|---|---|
| γ_φ' (tanφ' üzerine) | 1.00 | **1.25** |
| γ_c' (efektif kohezyon) | 1.00 | **1.25** |
| γ_cu (drenajsız kayma dayanımı) | 1.00 | **1.40** |
| γ_qu (serbest basınç dayanımı) | 1.00 | **1.40** |
| γ_γ (birim hacim ağırlık) | 1.00 | 1.00 |

### 1.3 Direnç (Resistance) — yüzeysel temel Tablo A.5 ✅
| Katsayı | R1 | R2 | R3 |
|---|---|---|---|
| γ_R;v (taşıma gücü) | 1.00 | **1.40** | 1.00 |
| γ_R;h (kayma) | 1.00 | **1.10** | 1.00 |

### 1.4 Direnç — genel/şev stabilitesi Tablo A.14 ✅
| Katsayı | R1 | R2 | R3 |
|---|---|---|---|
| γ_R;e | 1.00 | **1.10** | 1.00 |

### 1.5 Tasarım yaklaşımları (Design Approaches) ✅
| DA | Kombinasyon | Set birleşimi | FEM'e uyum |
|---|---|---|---|
| **DA1-C1** | Komb.1 | A1 + M1 + R1 | yapısal-yük baskın (malzeme azaltılmaz) |
| **DA1-C2** | Komb.2 | A2 + M2 + R1 | **malzeme-katsayılı → SRM** |
| **DA2** | — | A1 + M1 + R2 | direnç-katsayılı → rapor katmanı |
| **DA3** | — | (A1 yapısal / A2 geoteknik) + M2 + R3 | **malzeme-katsayılı → SRM** |

Kullanım: DA2 çoğu Avrupa ülkesinde taşıma gücü; DA3 çoğu ülkede şev/genel stabilite; DA1 UK/Portekiz.

## 2. TBDY 2018 Bölüm 16 — deprem etkisinde temel/zemin tasarımı

Kaynaklar: TBDY 2018 (Resmî Gazete, 18.03.2018); İMO Bölüm-16 sunumu; dergipark/akademik özetler.

### 2.1 Çerçeve ✅ (paradigma değişimi)
- **"Zemin emniyet gerilmesi / güvenlik sayısı" kavramı TERK EDİLDİ.** Yerine **tasarım dayanımı**:
  karakteristik dayanım `R_k`, bir **dayanım katsayısı** `γ_Rv`'ye bölünerek tasarım dayanımı
  `R_d = R_k / γ_Rv` elde edilir; kontrol **`E_d ≤ R_d`** (tasarım etkisi ≤ tasarım dayanımı).
  Bu, EC7 DA2 tipi **direnç-katsayılı** yaklaşımdır (LRFD benzeri).
- Bölüm 16 başlığı: *"Deprem Etkisi Altında Temel Zemini ve Temellerin Tasarımı İçin Özel Kurallar."*
- Statik **ve** deprem yük durumları ayrı ayrı kontrol edilir (`E_d` her iki durum için).
- **16.6 + Ek-16B:** basitleştirilmiş **sıvılaşma** tetiklenme potansiyeli değerlendirmesi.
- Eğimli arazide göçme riski varsa **şev stabilitesi** analizi zorunlu.

### 2.2 Dayanım katsayıları γ_Rv — Tablo 16.2 (Bölüm 16.8.2) ✅ (çapraz-doğrulandı; birincil teyit bekliyor)
> **DÜZELTME (2026-07-13):** γ_Rv, Tablo **16.1**'de DEĞİL (o site-sınıfı, §2.3) — **Tablo 16.2,
> Bölüm 16.8.2**'de. Değerler iki bağımsız kaynak + EC7 DA2 ile çapraz-doğrulandı; kesin kilit için
> birincil metin (kullanıcı PDF, Tablo 16.2) teyidi bekleniyor.
> Tasarım kontrolü: **`E_d ≤ R_d = R_k/γ_Rv`** (temel için `q_0 ≤ q_t = q_k/γ_Rv`). Aynı γ_Rv statik
> ve deprem için geçerlidir (fark yük birleşiminde + sismik taşıma-gücü formülünde, katsayıda değil).

| Mekanizma | Sembol | γ_Rv | EC7 DA2 (R2) karşılığı |
|---|---|---|---|
| Taşıma gücü (bearing) | γ_RV | **1.40** | γ_R;v = 1.40 ✓ |
| Kayma (sliding / friction) | γ_RH | **1.10** | γ_R;h = 1.10 ✓ |
| Pasif direnç (passive) | γ_RP | **1.40** | γ_R;e = 1.40 ✓ |

**ÖNEMLİ:** TBDY 2018 geoteknik direnç-katsayıları sayısal olarak **EC7 DA2 ile aynıdır** → B1.4 rapor
katmanı ikisini tek mekanizmayla (`R_d = R_k/γ_Rv`, `E_d ≤ R_d`) ele alabilir. Statik/deprem farkı
yalnızca yük birleşimi + sismik taşıma-gücü (atalet) düzeltmelerindedir.
Kaynaklar: TBDY 2018 Tablo 16.2, Bölüm 16.8.2 (istinatduvari.com; sanalsantiye.com); Resmî Gazete 18.03.2018.

### 2.3 Yerel zemin sınıfları — Tablo 16.1 ✅ (birincil TBDY 2018 metni, kullanıcı 2026-07-13)
Deprem tasarım spektrumu bu sınıfa göre türetilir (D/dinamik iş kolu). Üst 30 m ortalamaları:

| Sınıf | Zemin cinsi | V_S30 [m/s] | N60,30 [darbe/30cm] | c_u,30 [kPa] |
|---|---|---|---|---|
| **ZA** | Sağlam, sert kayalar | > 1500 | – | – |
| **ZB** | Az ayrışmış, orta sağlam kayalar | 760–1500 | – | – |
| **ZC** | Çok sıkı kum/çakıl/sert kil veya ayrışmış çatlaklı zayıf kaya | 360–760 | > 50 | > 250 |
| **ZD** | Orta sıkı–sıkı kum/çakıl veya çok katı kil | 180–360 | 15–50 | 70–250 |
| **ZE** | Gevşek kum/çakıl veya yumuşak–katı kil; ya da PI>20 & w>%40 ile toplamda 3 m'den kalın yumuşak kil (c_u<25 kPa) | < 180 | < 15 | < 70 |
| **ZF** | Sahaya-özel araştırma gerektiren: (1) sıvılaşabilir/yüksek-hassas kil/göçebilir zayıf-çimentolu; (2) >3 m turba ve/veya yüksek-organik kil; (3) >8 m yüksek plastisite (PI>50) kil; (4) >35 m çok kalın yumuşak/orta-katı kil | özel | özel | özel |

### 2.4 Deprem tasarım spektrumu ⏳ (D/dinamik iş kolu)
- **TDTH** (Türkiye Deprem Tehlike Haritaları) → `S_S` (kısa periyot) ve `S_1` (1.0 s) spektral
  ivme katsayıları; yerel zemin sınıfı (Tablo 16.1, §2.3) → `F_S`, `F_1` katsayıları → tasarım spektrumu.
- Tam spektrum şekli + F_S/F_1 katsayı tabloları D (dinamik) iş kolunda ayrı kilitlenecek.

## 3. v0.3 uygulama planı (B1)

- **B1.0** (bu belge): EC7 kilitlendi ✅; TBDY çerçeve ✅ + γ_Rv/spektrum ⏳ pinleme bekliyor.
- **B1.1:** kısmi-katsayı **veri modeli** (kod + set A/M/R + ulusal-ek; config-driven, EN default).
- **B1.2 (B1a):** EC7 DA1-C2 & DA3 → mevcut `strength_reduction.hpp` sabit-katsayı uzantısı
  (`c_d=c/γ_c`, `tanφ_d=tanφ/γ_φ`), yüklere A-set; SRM makinesi yeniden kullanılır.
- **B1.3 (B1b):** direnç-katsayılı kontrol katmanı (EC7 DA2 + TBDY `E_d ≤ R_d`) — FEM çıktısından
  mekanizma direnci → `/γ_R` → geçti/kaldı.
- **B1.4:** kod-referanslı rapor satırı (hangi katsayı hangi maddeden) + V&V testi
  (katsayı=1.0 → mevcut analizle birebir; bilinen DA-set → kapalı-form kontrol).

## 4. Doğrulama (DoD, roadmap §6)
- γ=1.0 dejenerasyonu mevcut analizle **bit-birebir** (regresyon güvencesi).
- DA1-C2 malzeme azaltması SRM ile tutarlı (bilinen c/φ azaltma = FoS makinesi).
- Kapalı-form kontrol: şerit temel taşıma gücü, azaltılmış c/φ ile Rankine/Prandtl analitiği.
