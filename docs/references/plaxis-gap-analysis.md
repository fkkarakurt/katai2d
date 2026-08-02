# PLAXIS 2D 2025.1 — Gap Analizi ve Rota (KATAI vs PLAXIS)

> ⚠️ **TARİHSEL BELGE (2026-06 durumu).** Buradaki birçok ❌ artık ✅'dır (embedded beam,
> Biot konsolidasyon, transient/kuplajlı akış, Dynamics lineer+NONLİNEER, HSsmall, nokta yük,
> fazlar-arası yapısal durum taşıma...). **Güncel kapsam matrisi ve rota: `docs/internal/ROADMAP.md` §3-5**;
> güncel doğrulama durumu: `docs/validation/VALIDATION-SUMMARY.md`. Bu dosya, PLAXIS manual
> bölüm-haritası referansı olarak korunuyor (bölüm numaraları hâlâ geçerli).

Kaynak: **PLAXIS 2D 2025.1 Scientific Manual** + **Reference Manual** (Seequent/Bentley, 24 Eyl 2025;
files.seequent.com/PLAXIS/Manuals/PLAXIS_2D/English/). İki manual'in TAM içindekiler tablosu çekildi/okundu.
Amaç: KATAI'nin neyi karşıladığını, neyin eksik olduğunu PLAXIS'in kendi yapısına göre belirlemek ve
önceliklendirmek. ✅=var+doğrulu, ◑=kısmi, ❌=yok. **Mimari teyidi:** Scientific Man. §2 deformasyon teorisi
(LᵀΣ+b=0, ε=Bv, Galerkin zayıf form, σⁱ=σⁱ⁻¹+Δσ) KATAI ile BİREBİR → temel doğru.

## A. Çözücü / deformasyon teorisi (Sci. Man. §2, §9, App B/C)
| Konu | PLAXIS | KATAI | Not |
|---|---|---|---|
| Continuum + FE discretisation (§2.1-2.2) | ✅ | ✅ | birebir aynı |
| Implicit plasticity integration (§2.3) | ✅ | ✅ | MC kapalı-form, HS substepping |
| **Global iterative method (App B/C)** | **elastik-rijitlik + Quasi-Newton + over-relaxation** | **tam-Newton (tutarlı teğet) + line search + adaptif adım** | ⚠ FARKLI ALGORİTMA. KATAI daha modern; PLAXIS yakınsama-yolunu birebir tutturmak istersek elastik-rijitlik seçeneği eklenebilir |
| Convergence criteria (§9) | global+local error (out-of-balance/CSP) belirli formül | residual-norm/ref | ◑ PLAXIS formülüne hizalama = "aynı kriter" |
| **Arc-length control (Ref §7.8)** | ✅ | ❌ | snap-through/limit yük için |
| Otomatik adım boyutu (Ref §7.8.1) | ✅ | ✅ | KATAI adaptif cutback |

## B. Elemanlar (Sci. Man. §7)
| Eleman | PLAXIS | KATAI |
|---|---|---|
| 6- ve 15-düğümlü üçgen | ✅ | ✅ tri6+tri15 |
| Plate (3/5-düğüm line) | ✅ | ✅ (tri6+tri15 kenarı) |
| Geogrid | ✅ | ✅ |
| Anchor (fixed-end + node-to-node) | ✅ | ✅ elastoplastik |
| Interface (Goodman) | ✅ | ✅ (+K0 seeding) |
| **Embedded beam (§7.5, pile row)** | ✅ skin+foot etkileşim | ❌ | formülasyon manual'de |
| Cables, discontinuities, connections, tunnels | ✅ | ❌ | ileri/özel |

## C. Malzeme modelleri (Material Models Man. + Ref §6)
| Model | PLAXIS | KATAI | Öncelik |
|---|---|---|---|
| Linear Elastic | ✅ | ✅ | — |
| Mohr-Coulomb | ✅ | ✅ doğrulu | — |
| Hardening Soil | ✅ | ✅ büyük oranda | — |
| **HSsmall** | ✅ | ✅ (yasa+FE, commit 1e4368f/49dc60f) | TAM — round-off vs MMM §7 |
| Soft Soil / Soft Soil Creep | ✅ | ❌ | orta |
| Modified Cam-Clay | ✅ | ❌ | orta |
| Hoek-Brown (kaya) | ✅ | ❌ | düşük |
| Undrained A / B / C (Ref §6.2) | ✅ | A✅ B✅ C❌ | B yeni doğrulandı |
| Unsaturated soil | ✅ | ❌ | düşük |

## D. Yeraltısuyu akışı (Sci. Man. §3)
| Konu | PLAXIS | KATAI |
|---|---|---|
| Steady-state confined/unconfined | ✅ | ✅ (Darcy/Laplace/uplift/freatik/sızma yüzü) |
| **Transient (zaman-bağımlı) akış (§3.1.1)** | ✅ | ❌ |
| Flow in interface elements (§3.4) | ✅ | ❌ |
| Wells, drains (Ref §5.9) | ✅ | ◑ (Neumann akı var) |

## E. Konsolidasyon (Sci. Man. §4 — Biot)
| Konu | PLAXIS | KATAI | Öncelik |
|---|---|---|---|
| **Biot kuplajlı konsolidasyon** | ✅ | ❌ (undrained anlık var, zaman-bağımlı yok) | **YÜKSEK (Terzaghi 1D klasik kontrol)** |
| Fully coupled flow-deformation (Ref §7.4.4) | ✅ | ❌ | orta |

## F. Analiz tipleri (Ref §7.4)
Initial stress/K0 ✅ · Plastic ✅ · **Consolidation ❌** · **Fully coupled ❌** · Safety (φ/c reduction) ✅ ·
Dynamics ❌ · Plastic nil-phase ◑ (gravity dengesi) · **Updated mesh (büyük deformasyon) ❌**

## G. Yükler ve BC (Ref §5.3-5.4)
Dağılı yük ✅ (surface traction) · **Nokta yük ◑ (doğrula/ekle)** · Prescribed displacement ✅ ·
Line contraction (tünel) ❌

## H. Staged construction (Ref §7.5, §7.11)
Eleman aktif/pasif (kazı/dolgu) ✅ · Çok-faz gerilme aktarımı ✅ · **Interface durum faz-taşıma ❌** ·
Anchor prestress ◑ · ΣMstage ✅ (staged-release ramp) · Volumetric strain ◑

## I. Soil lab test simülasyonu (Ref §11.12)
Triaxial/oedometer KATAI'de malzeme-testi olarak ✅ (HS kalibrasyon); kullanıcıya "SoilTest" aracı ❌ (GUI).

## J. Düşük öncelik (ileri/özel)
Thermal/THM (Sci §5), Dynamics (Sci §6), Sensitivity/parameter variation (Sci §8), Design approaches (Ref §5.8),
tüneller, GUI/Output programı, Paraview export.

---

## ROTA — önceliklendirilmiş (insanların ilk bakacağı + çekirdek + manual referanslı)
**Faz A (en yüksek görünürlük, çekirdek doğruluk):**
1. **HSsmall** — Material Models Man.; HS üstüne G₀/γ₀.₇ (Hardin-Drnevich). Yeni mimari reçetesinin ilk testi.
2. **Biot konsolidasyon** — Sci. Man. §4; Terzaghi 1D U-Tv analitik doğrulama (altın standart <%2). Ek u-DOF (pore).
3. **Nokta yük** — Ref §5.3.1; Flamant analitik ile doğrula (kolay, eksik kapat).
4. **Embedded beam** — Sci. Man. §7.5 (formülasyon hazır); pile row skin+foot.

**Faz B (çözücü/algoritma PLAXIS-uyumu):**
5. **Convergence criteria** PLAXIS §9 formülüne hizala + **arc-length** (Ref §7.8) — "aynı algoritma".
6. **Transient groundwater flow** (Sci §3.1.1) + flow in interfaces.
7. Wall benchmark: **graded/unstructured mesh** (Ruppert) → derin vakalar <%5 (domain verimli).

**Faz C (kapsam genişleme):**
8. Soft Soil / Modified Cam-Clay; Undrained C; fully coupled flow-deformation; updated mesh; dynamics.

**Doğrulama disiplini her adımda:** önce matematik referans (PLAXIS manual + makale), izole test → analitik
kapalı-form (<%1-2) → kuplajlı → çok-kod bandı ortası. Kurulum doğruluğu (mesh/domain/BC/faz) = en sık hata.

İlgili: [[material-model-architecture]], [[multi-code-validation-plan]], [[validation-benchmarks]], [[literature-review]].
