# GUI compute-path (entegre `build_problem`) doğruluk denetimi

Bu doküman, KATAI'nin **GUI'nin gerçekte kullandığı uçtan-uca yolu** (`studio/app/build_problem.hpp` →
`mesh_from_project` → `solve_gravity_le` → nodal stress recovery) kapalı-form/PLAXIS referanslarıyla
**niceliksel** doğrular. `VALIDATION-SUMMARY.md` çekirdek kernel'leri *doğrudan* doğrular (yapısal mesh,
elle kurulan DOF'lar); buradaki amaç farklı: **kullanıcının GUI'de bastığı "Calculate"** aynı doğru sonucu
veriyor mu? "Pre + compute + post"un tamamı denetlenir.

V&V disiplini (Roache 1998): en basitten karmaşığa. Önce lineer-elastik (kapalı-form tam), sonra
nonlineer (Mohr-Coulomb taşıma gücü), sonra su/yapısal/staged.

## 0. Düzeltilen kusur: K0 procedure POST-PROCESS gösterimi (2026-06-10, commit 74057ec)
Kullanıcı "GUI'de K0 procedure çalışmıyor gibi" bildirdi. **`diag_k0` ile kanıtlandı: K0 procedure DOĞRU
HESAPLIYOR** — sorun yalnız gösterimdeydi:
- **K0 max|u| = 0.0 (TAM)**, `σ'_v = −γ(H−y)` analitikle **1e-13**, `σ'_h = K0·σ'_v` tam.
- Gravity loading 3.34 cm oturur ama **aynı** `σ'_v`'ye varır (düşey denge yoldan bağımsız) — ikisi yalnız
  deplasmanda farklı (PLAXIS'in tam ayrımı).

İki gösterim hatası düzeltildi: (1) varsayılan alan `u_y` idi → K0'da ~0 → "boş/çalışmıyor" görünüyordu;
artık K0 solve sonrası varsayılan **`σ'_yy`** (PLAXIS de K0 sonrası efektif gerilme gösterir). (2) deplasman
auto-ölçeği `max|u|`'ye bölüyordu → K0'ın ihmal-edilebilir deplasmanını round-off patlaması yapıyordu;
artık `|u| < 1e-7·H` ⇒ deforme gösterilmez.

## 1. Lineer elastik (kapalı-form) — GUI yolundan TAM
| Problem | KATAI (GUI) | Referans | Hata | Test |
|---|---|---|---|---|
| 1B ödometre öz-ağırlık `u_top=−γH²/2E_oed` | FE | kapalı-form | <%2 | `test_solve` |
| 1B ödometre dağılı surcharge `−qH/E_oed` | FE | kapalı-form | <%2 | `test_solve` |
| K0 procedure max\|u\| | 0.0 (tam) | 0 (jeostatik denge) | round-off | `test_gui_solve`, `diag_k0` |
| K0 efektif `σ'_v=−γ(H−y)`, `σ'_h=K0σ'_v` | FE | kapalı-form | 1e-13 | `test_gui_solve`, `diag_k0` |
| Batık blok buoyant `σ'_v=−γ'(H−y)` + pore `γ_w(H−y)` | FE | Terzaghi | round-off | `test_gui_solve` |

## 2. Mohr-Coulomb — GUI yolundan niceliksel (`study_gui_validation`, on-demand ~45s)
| Problem | KATAI (GUI) | Referans | Hata | Not |
|---|---|---|---|---|
| Prandtl şerit temel `N_c` (φ=0, ağırlıksız) | 4.84 | 2+π = 5.14 | **−%5.9** | artımsal limit-analiz (alttan/konservatif) |
| Reissner `N_c` (φ=20°, ağırlıksız) | 13.94 | 14.84 | **−%6.1** | (N_q−1)cot φ |
| MC = LE elastik rejimde (ödometre, akma-altı) | birebir | LE `−qH/E_oed` | <0.1% (MC vs LE) | sahte akma/yumuşama YOK |

**Çöküş-yükü bulgusu (kusur değil, karakteristik):** GUI yolu çöküş yükünü ~%6 **alttan** (konservatif,
güvenli taraf) tahmin ediyor. Sebep artımsal limit-analizin **yük-adımı granülaritesi**: GUI'nin nonlineer
şeması 20 yük adımı kullanır (`min_dlam ≈ 1/160`), çöküşe tam yaklaşamadan en yüksek dengeli seviyede durur.
İthaf edilmiş çekirdek yol (`test_prandtl`, 40 adım + yapısal mesh) ~%1 alır. **Elastik/servis tepkisi TAM.**
Daha ince mesh çöküşü İYİLEŞTİRMEZ (mekanizma daha erken oluşur, adım granülaritesi limit kalır) — bu yüzden
çöküş-yükü doğruluğu mesh değil adım-çözünürlüğü ile sınırlı. (PLAXIS de SumMstage ile çöküşe yaklaşır.)

`SolveResult.load_factor` eklendi (commit ek): dengelenen yük oranı = artımsal limit yük (PLAXIS SumMstage);
GUI artık yakınsamamada "yükün %X'i dengelendi (çöküş mekanizması)" der, ham "did not converge" yerine.

## 3. Axisymmetric (PLAXIS 2. mod) — GUI yolundan TAM (commit 56cbdb4)
Çekirdek axisym tamamen doğrulanmıştı (`test_axisym_*`, Lamé 1e-5, silindir çökme %0.3) ama **GUI'den
erişilemiyordu** — denetimin ortaya çıkardığı en büyük "güç var ama bağlı değil" boşluğu. Artık
`project.axisymmetric` ile bağlı: r-ağırlıklı gravity/internal-force/traction + `Kinematics::Axisymmetric`.
v1 kapsamı zemin-only + kuru (su/yapısal axisym = dürüst guard mesajı). `test_axisym_gui` (hepsi TAM):
| Problem | KATAI (GUI axisym) | Referans | Hata |
|---|---|---|---|
| K0 silindir (max\|u\|, σ_zz, σ_rr) | 0.0 / FE | jeostatik denge | 1e-13 |
| Gravity ödometre `−γH²/2E_oed` | −0.066859 | −0.066857 | <0.01% |
| Tek-eksenli oturma `qH/E` | −0.020000 | −0.020000 | round-off |
| **Radyal genişleme `u_r(R)=νqR/E`** | **0.003000** | **0.003000** | **round-off** |

Son satır KRİTİK: hoop şekildeğiştirmesi `ε_θ=u_r/r` axisym'i düzlem-şekildeğiştirmeden AYIRAN özellik;
birebir doğru → r-ağırlıklı montaj + 4-bileşenli bünye + hoop kuplajı entegre yolda DOĞRU.

## 4. Sistematik matris — structures × yük × su × model × tri6/tri15 (`test_gui_matrix`)
Kullanıcı direktifi: "structures, yük ve yer altı suyunun hesaba katıldığını, her zemin modeli, tri6/tri15
ve her yapı elemanı için kontrol et." Referans-bağımsız ama KESKİN değişmezler, tam GUI yolundan (~3s):
- **(I1) K0 admissibility — 20/20 TAM:** LE+MC × tri6+tri15 × {none, plate, anchor, geogrid, embedded
  beam} → öz-ağırlık K0'da `max|u|=0` (tam). En keskin pre+compute kontrolü: yanlış K0/yapı-montaj/pre
  etkileşimi sahte "kurulum" hareketi üretir. Hepsi sıfır → K0 + tüm yapılar + iki düzen + LE/MC TUTARLI.
- **(I2) Her yapı hesaba katılıyor:** 300 kN yük altında plate %80.8, anchor %18.4 (fixed-end), geogrid
  %3.0, embedded beam %11.9 değiştirir (montaj ediliyor) + alan sonlu/sınırlı.
- **(I3) tri6 ≈ tri15:** öz-ağırlık (tekil-olmayan) oturması **%0.00** (birebir) → tri15 GUI yolu doğru.
  (NOT: nokta yükü TEKİLdir → settlement ıraksar, derece-bağımlı; tri6/tri15 kıyası tekil-olmayan alanla.)
- **(I4) Su buoyancy:** tri6+tri15 taban `σ'_v=−γ_unsat·2−γ'·8` (kısmi su tablası) <%0.05.

## 5. Safety analizi (φ-c reduction / SRM) — GUI yolundan TAM (commit sonraki)
Kullanıcı: "K0'da |u|=0 imkansız, zeminin ağırlığı var; slope'ta çökme görmem gerekmez mi?" → KAVRAMSAL
açıklama + EKSİK ÖZELLİK. **K0'da |u|=0 DOĞRU** (PLAXIS gibi): K0 procedure başlangıç GERİLMESİNİ kurar,
deplasmanı değil — bozulmamış zemin jeolojik olarak zaten kendi ağırlığı altında dengede (gerilme
σ'_v=γz ağırlığı YANSITIR; deplasman bu dengeye GÖRE sıfır). Çökmeyi görmek için ayrı analiz gerekir:
**Safety (φ-c reduction)** — çekirdek `factor_of_safety` doğrulanmıştı (test_slope FoS 1.01) ama GUI'de YOKTU.
Bağlandı: `safety_analysis` (çok-malzeme c/φ azaltma + mekanizma NewtonResult) + build_problem `InitialPhase::
Safety` + GUI "Safety (phi-c reduction)" + FoS gösterimi. `test_safety_gui` (Griffiths&Lane 1:2 şev):
**FoS(GUI)=1.009 vs ~0.99 (%1.9)**, mekanizma max|u|=1.13 (slip yüzeyi GÖRÜNÜR), LE dürüstçe reddedilir.
Üç analiz: K0=başlangıç gerilmesi (u≈0), Gravity loading=öz-ağırlık oturması, Safety=çökme/FoS.

## 6. K0 eğimli yüzeyde — denge (nil) adımı (2026-06-11, commit 8c0e421) — `test_k0_slope`
PLAXIS kuralı: K0 yalnız yüzey+katman+su tablası YATAYKEN doğru; değilse gerçek dengesiz kuvvetler →
plastik nil-step. KATAI: dengesizlik `d = f_body − f_int(seed)` dış yüklerle ramplanır (tetik GEOMETRİK,
|d| eşiği değil — kuadratür artığı nil-step taklidi yapamaz). Ayrıca: kolon integrali strata-break
KESİN (katmanlı/su'lu yatayda u=0 korunur) + sınır düğümünde çakışık kenarların **BC birleşimi**
(Free kenar sabit kenarı ezemez — katman arayüzü × yan sınır mesnet deliği kapatıldı).
| Kontrol | Sonuç |
|---|---|
| Yatay 1-katman / 2-katman / yatay su: max\|u\| | **0.0 TAM**, nil-step yok |
| Griffiths&Lane şev (stabil c=20): nil-step | aktif; max\|σ_xy\| ≈ 38 kPa gelişir |
| Kret altı σ'_v vs overburden / vs gravity loading | < %10 / < %10 (PLAXIS "karşılaştırılabilir") |
| Stabil-olmayan şev (c=1) | dürüst çöküş: load_factor 0.83, ok=false |
| Eğimli su tablası | nil-step tetiklenir |

## 7. Groundwater flow GUI yolu (2026-06-11) — `test_seepage_gui`
`build_flow.hpp`: kenar-başına akış BC (Closed/Head/Seepage face, sağ-tık) → konfine HIZLI-YOL (tek
lineer çözüm; ψ≥0 her yerde ise k_rel≡1 zaten kesin) ya da unconfined/seepage-face aktif-set Picard.
Mekanik kuplaj: `solve_gravity_le(..., flow_head)` — pore + doygunluk AKIŞ head alanından
(`assemble_pore_load_from_head` + YENİ `assemble_gravity_from_head`, aynı interpolasyon = tutarlı çift);
K0 fazı akış kuplajında daima nil-step (seed hidrostatik kalır, hedef gerçek denge).
| Problem | KATAI (GUI) | Referans | Hata |
|---|---|---|---|
| Konfine Darcy 1B (iki rezervuar) max\|h−h_ex\| | 1.6e-13 | lineer kapalı-form | **round-off** |
| Konfine debi Q = k·i·A | 1.500000 | 1.5 m³/gün/m | **kesin** |
| Charny seepage-face debisi q=k(h₁²−h₂²)/2L | 0.612 | 0.600 | %2.0 |
| Terzaghi yukarı-sızma σ'_v taban (kuplajlı K0) | −82.280 | −82.280 | **<1e-3** |
| Pore gösterimi (akış head'inden) | 117.720 | 117.720 | kesin |
| Kütle dengesi Σ Q | ~1e-13 | 0 | round-off |

## 8. Yapısal M/Q/N çıktıları GUI'de (2026-06-11) — `test_struct_gui`
`SolveResult.struct_forces`: çizilen her yapısal hat için diyagram (structural_forces.hpp
post-processor'ları — çekirdek matematik test_struct_forces'ta kapalı-formla 2.7e-9'da kilitli;
buradaki denetim WIRING + fizik). GUI: Results panelinde zarf + M/Q/N eğri (PlotLines) + tablo;
viewport'ta PLAXIS-tarzı eksene-dik ofset diyagram overlay'i (M kırmızı / Q yeşil / N mavi, tepe
değer etiketli); anchor N (kapasite uyarısıyla).
| Kontrol | Sonuç |
|---|---|
| Yatay radye, merkez yükü: M piki | s = 4.00 / L = 8.00 (tam orta-açıklık) |
| Serbest radye uçlarında M | zarfın %2.7'si (≈0, moment-free uç) |
| Fixed-end tie 100 kN yük altında | N = 47.5 kN/m çekme (0 < N < P payı) |
| Geogrid | tension-only N ≥ 0, istasyonlar hat üstünde |
| **Embedded beam (pile row) N/Q/M (2026-06-13)** | `embedded_beam_force_diagram` (plate::forces kernel'i yeniden — beam'in node_x/y + dof_x/y/phi); kafa eksenel yükte **N(head)=P %0.07**, serbest ayağa sönüm (|N| ayak %0.3), M≈0; GUI-path sürşarjla pile hattı (x=10) max\|N\|=260 |
| İstasyon konumları / yay-uzunluğu | çizilen hat üstünde, monoton |

## 9. Yerel mesh yoğunluğu (2026-06-11) — `test_mesh_refine`
Sizing-field'li Ruppert + PLAXIS coarseness semantiği (bölge/hat/yük f, refine ×0.5 / coarsen ×2,
EMR-vari auto-refine GUI'de açık) — formülasyon + tablo: `docs/references/mesh-sizing.md`.
Factor'sız model eski mesher ile bit-birebir; bölge f=0.25 alan cap'i tam; hat kaynağı gradasyonlu;
min-açı 22° korunur. GUI: Mesh panelinde min-açı/ortalama-kalite istatistikleri, sağ-tık refine menüsü.

## 10. Çok-fazlı staged construction (2026-06-11) — `test_staged_gui`
Faz yöneticisi (PLAXIS Phases): her faz konfigürasyon dengesizliğini ramplar —
`ramp = f(aktif) − f_int(committed)` = **SumMstage** (kazı boşaltması, dolgu ağırlığı, yeni yük
TEK kuraldan). Initial faz aktivasyonu (kısmi geometri K0), committed Gauss zinciri, pasif eleman
gerilme-kurtarma dışlama, Safety fazı (aktif konfig FoS). GUI: faz listesi + ekle/sil (önceki fazın
konfigini miras alır), faz başına aktif nesne checkbox'ları, faz sonucu seçici + kümülatif deplasman.
| Problem (1B kapalı-form, E_oed) | KATAI (GUI) | Kesin | Hata |
|---|---|---|---|
| Initial (yalnız alt tabaka aktif, K0) max\|u\| | 0.0 | 0 | **kesin** |
| Dolgu: eski yüzey oturması γ_f·h_f·H/E_oed | −0.03031 | −0.03031 | **<1e-4** |
| Dolgu üstü (öz-sıkışma dahil) | −0.04041 | −0.04041 | **<1e-4** |
| Dolgu sonrası taban σ'_v | −176.00 | −176.00 | **kesin** |
| 3. faz sürşarj artımı q·H/E_oed | −0.03714 | −0.03714 | **<1e-4** |
| Kazı: taban kabarması (boşaltma) | +0.03031 | +0.03031 | **<1e-4** |
| Kazı sonrası taban σ'_v | −108.00 | −108.00 | **kesin** |

## 11. Konsolidasyon (zaman-bağımlı, Biot) GUI yolu (2026-06-13) — `test_consolidation_gui`
PLAXIS "Consolidation" hesap fazı (`PhaseType::Consolidation`): fazın yük artımı dF = f − B (SumMstage)
t=0+'da uygulanır → S≈0 (gerçek su Kw/n) ⇒ UNDRAINED tepki fazlalık boşluk basıncı üretir, sonra `duration`
[gün] boyunca sönümlenir (oturma gelişir). Drenaj sınırı Flow conditions FlowBC'den (Head/Seepage = açık);
`active` maskesi çekirdeğe geçer (kazı/dolgu sonrası konsolidasyon). Çözücü **sparse PARDISO sym-indefinite
(mtype=−2) saddle-point**, factor-once-solve-many (dense=sparse birebir doğrulandı). Çıktı: oturma-zaman +
fazlalık-pore-zaman eğrisi + nihai alanlar.
Yanal-konfine kolon, üstte drenajlı, üst sürşarj q: s∞ = qH/E_oed, U(Tv) = Terzaghi serisi, Tv = c_v t/H², c_v = k E_oed/γ_w.
| Problem (1B Terzaghi, GUI yolu) | KATAI (GUI) | Kapalı-form | Hata |
|---|---|---|---|
| Undrained fazlalık pore üretimi (t=0+, Skempton B) | 9.99 kPa | q = 10 | **%0.1** |
| Nihai fazlalık pore (Tv≈2) | 0.10 kPa | →0 | **<%1** |
| U_FE(Tv=0.2) | 0.4982 | 0.5041 | **−%1.2** |
| U_FE(Tv=0.6) | 0.8100 | 0.8156 | **−%0.7** |
| U_FE(Tv=0.9) | 0.9080 | 0.9120 | **−%0.4** |
| Nihai oturma (Tv≈2) U·s∞ | 0.11923 m | 0.11930 | **%0.06** |

## 12. Undrained (A/B) drenaj iş akışı GUI yolu (2026-06-13) — `test_undrained_gui` + `study_gui_validation`
PLAXIS drenaj tipi (Material > General). **Undrained (A):** efektif dayanım c',φ' + pore-akışkan rijitliği Kw/n
(νu=0.495) solver teğetine + toplam gerilmeye girer (efektif gerilme yolu; fazlalık pore eps_vol'de). **Undrained
(B):** aynı Kw/n makinesi ama undrained dayanım doğrudan girilir — c = su, φ ZORLA 0 (Tresca su; build_problem MC
için φ=ψ=0 yapar). Enum dosya-uyumlu (UndrainedB=3 eklendi); GUI PLAXIS sırasında index↔enum haritasıyla sunar.
Konfine kolon (ağırlıksız, üst sürşarj q): M_u=M'+Kw/n, oturma −qH/M_u, efektif σ'_yy=−M'q/M_u, kalan ~q fazlalık pore.
| Problem (GUI yolu) | KATAI (GUI) | Kapalı-form | Hata |
|---|---|---|---|
| Undrained (A) konfine oturma | −qH/M_u | −qH/M_u | **<%2 (birebir)** |
| Drained konfine oturma | −qH/M' | −qH/M' | **<%2 (birebir)** |
| M_u/M' (νu=0.495, ν'=0.3) → fazlalık pore u/q | 28.9× / 0.965 | Skempton ~1 | **doğru** |
| **Undrained (B) şerit temel q_ult (φ'=30 GİRİLİR, yok sayılır)** | **5.175 su** | **(2+π) su = 5.142** | **+%0.6** |
*(Undrained B doğrulaması KRİTİK: girilen efektif φ'=30 ZORLA 0'a iner → 5.14 su; sızsaydı ~5× büyük çıkardı.
study_gui_validation on-demand limit analizi; bearing_q_ult drainage parametresi.)*

## 13. Kalan denetim (sıradaki, foundations-up sırası)
- **HS GUI yolundan niceliksel** (tek-eleman triaxial/oedometer GUI'den; yüklü BVP robustluğu — bilinen iş).
- ~~Konsolidasyon elastoplastik (MC/HS)~~ **TAM (2026-06-13)** — monolitik kuplajlı-Newton (test_consolidation_plastic: LE-indirgeme 1e-16, MC footing drenajlı-limit %0.3); build_problem MC/HS'de plastik yolu çağırır. KALAN: HS undrained-B su-cap konsolidasyon kalibrasyonu.
- **Undrained (C) total-stress** (efektif gerilme yok, νu+su; v1 yok — A/B lokmasız, C lokmaya yatkın) + HS undrained-B (su cap).
- **Axisym genişletme:** su (r-ağırlıklı pore) + yapısal (shell/ring) — v2; axisym + flow kuplajı da guard'lı.
- **Post-process gerilme geri-kazanımı doğruluğu** (SPR; lineer alanda 1e-13 kanıtlı, eğri alanda kontrol).
- **Flow + diyagram + faz GUI etkileşimi göz-doğrulaması** (BC işaretleri, head/pore render, M/Q/N overlay,
  faz listesi/aktivasyon) — compute yolu test edildi, görsel yol manuel kontrol bekliyor.
- ~~Embedded beam (pile row) kuvvet diyagramı~~ **TAM (2026-06-13)** — embedded_beam_force_diagram (plate::forces yeniden); build_problem diag kind=3 → SolveResult.struct_forces; GUI N/Q/M overlay/tablo (kind 0). test_embedded_beam + test_struct_gui.
- **Safety fazında yapısal kuvvetler** — mekanizma anındaki (azaltılmış dayanım) kuvvetler raporlanmıyor (bilinçli).
- **Faz v1 sınırları (dürüst guard'lı):** gömülü perde (split) faz başına aktif/pasif edilemez; akış kuplajı
  fazlar içinde yok; yapısal plastik durum (anchor U_p, geogrid ε_p) fazlar arası taşınmıyor; axisym fazsız.
- **Kazılmış bölge render'ı** — pasif elemanlar alanda sıfır görünür (gri/gizli çizim katai_render işi).

## 14. SAĞLAMLIK (rock-solid) sözleşmesi — `test_robustness` (2026-06-13)
Ticari ürün kuralı: son kullanıcının KARMAŞIK veya bozuk modelinde program **ÇÖKMEMELİ**. GUI hesap
giriş noktaları (`mesh_from_project`, `solve_gravity_le`, `solve_phases`, `solve_groundwater_flow`)
GUI tarafından try/catch OLMADAN çağrılır → sözleşme: **asla fırlatma, asla hang, asla NaN/Inf döndür** —
daima `{ok=false + dürüst mesaj}` ya da sonlu çözüm. `test_robustness` 51 saldırgan/dejenere/aşırı
konfigürasyonu tüm hatta fırlatır (boş/atanmamış malzeme; sıfır-alan/collinear/bowtie/dup/NaN geometri;
E=0/<0/NaN/Inf/1e308, ν=0.5/0.7/−3; BC'siz rijit cisim; Inf/NaN yük; undrained ν'=0.5 Kw/n tekil;
konsolidasyon duration/steps=0/<0, kx=0; faz tüm soil pasif; flow BC'siz/k=0; **yapısal: sıfır-uzunluk/
NaN/domain-dışı plate/embedded-wall/anchor/geogrid/beam, yatay perde**; axisym negatif yarıçap; +
**KARMAŞIK geçerli model: 2-katman+gömülü perde+strut+sürşarj+su+kazı fazı**). **SONUÇ: 51/51 hayatta
kaldı — sıfır çökme, sıfır NaN** (düzeltme gerekmedi; önceki oturumların guard'ları zaten kapsamlı).
Değer: bu davranış artık REGRESYON olarak kilitli — yeni özellik bu girdilerde çökme sokarsa test yakalar.
Bilinen yumuşak nokta (çökme DEĞİL, gelecek input-validation): bazı dejenere konfig (bowtie, φ=120°,
axisym x<0) reddedilmek yerine sonlu-ama-anlamsız çözüyor.

Çalıştırma: `cmake --build build/msvc-rwdi --target study_gui_validation && bin/study_gui_validation`
(taşıma gücü); `bin/test_robustness` (çökme-dayanıklılığı); `diag_k0` (K0 denetimi). Kaynaklar: Prandtl 1921,
Reissner 1924, Terzaghi 1943, Roache 1998 (V&V), PLAXIS 2D Reference Manual (K0 procedure, SumMstage).
