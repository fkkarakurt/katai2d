# Malzeme Modeli Mimarisi — Yeni Zemin Modeli Eklemek İçin Sistem

Hedef: yeni bünye (constitutive) modellerini **kolay, güvenli ve sektör programlarıyla uyumlu**
şekilde eklemek. KATAI'nin malzeme katmanı, PLAXIS 2D / GEO5 / FLAC / Midas GTS NX ile **aynı
sayısal mimariyi** kullanır: **gerilme-noktası (Gauss) seviyesinde return-mapping + tutarlı
(algoritmik) teğet**. Bu doküman mimariyi ve "yeni model nasıl eklenir" reçetesini kilitler.

## 1. Çekirdek arayüz — malzeme noktası (material point)
Tüm bünye modelleri tek bir imza üzerinden çağrılır (`materials/material_model.hpp`):
```
integrate_point(const MaterialModel& m, const GaussState& committed,
                const StrainIncrement& dε, GaussState& trial, Tangent& D_T)
```
Sözleşme (PLAXIS/FLAC/Abaqus UMAT ile birebir aynı kavram):
- **Girdi:** yakınsamış durum `committed` (σ_n + iç değişkenler) ve şekildeğiştirme artımı `dε`.
- **Çıktı:** `trial` (σ_{n+1} + güncel iç değişkenler) ve **tutarlı teğet** D_T = ∂σ_{n+1}/∂dε.
- **Saf fonksiyon:** `committed`'in saf fonksiyonu (yol-bağımlı durum committed'te) → line-search ve
  adaptif yük adımı güvenli (ara denemeler committed'i bozmaz).

Bu, klasik **elastik öngörücü → plastik düzeltici (return mapping)** kalıbıdır (Simo & Hughes 1998;
Potts & Zdravković 1999). Çözücü (`nonlinear_solver`) modelden BAĞIMSIZDIR — yalnız bu arayüzü çağırır.

## 2. Durum (state) ve parametreler
- **`GaussState`** (yol-bağımlı iç değişkenler): `stress` (σ_xx,σ_yy,σ_xy), `stress_zz`, `eps_vol`
  (undrained excess pore), `gamma_p` (HS shear hardening), `pp` (HS cap preconsolidation). Yeni model
  yeni iç değişken gerektirirse buraya eklenir (committed/trial solver tarafından taşınır).
- **`MaterialModel`** (sabit parametreler, etiketli/flat — vtable YOK, hot-loop için): `type`
  (`MaterialType` enum), E, ν, c, φ, ψ, undrained bayrağı, `HardeningSoilParams hs`, ... Yeni model
  parametreleri ya doğrudan ya da bir alt-struct (hs gibi) olarak eklenir.

## 3. Kinematik bağımsızlık (plane strain + axisymmetric + ...)
Return mapping **asal-gerilme uzayında** yazılır → kinematik (plane strain / axisymmetric) bağımsızdır.
Çözücü `Kinematics` politikası ile templatelenir (`PlaneStrainKin`/`AxisymKin`): B matrisi, gerilme
boyutu, pore-vektörü politikadan gelir; return mapping aynı kalır. (MC bunu kanıtladı: aynı
`mc_return_mapping` hem plane strain hem axisym'de kullanılıyor.) Yeni model otomatik her iki modda çalışır.

## 4. Yeni bir zemin modeli eklemek — reçete
1. **Matematiği kilitle** `docs/references/<model>-formulation.md` (akma yüzeyi, akış kuralı, sertleşme,
   teğet; KAYNAK: PLAXIS MMM / FLAC / GEO5 / orijinal makale). KATAI kuralı: önce referans, sonra kod.
2. **Çekirdek model dosyası** `materials/<model>.hpp` — saf, basınç-veya-tension konvansiyonu net;
   `<model>_return_mapping(predictor, params) -> {stress, plastic, tangent, new_state}`. İzole test edilebilir.
3. **Parametre + durum:** `MaterialType::<Model>` enum'a; parametreleri `MaterialModel`'e (veya alt-struct);
   yeni iç değişken varsa `GaussState`'e.
4. **Dispatch:** `integrate_point` (plane strain) ve `integrate_point_axisym` switch'lerine `case` ekle
   (asal-uzay return mapping ikisinde de yeniden kullanılır → genelde tek helper).
5. **Tutarlı teğet:** kapalı-form (tercih) veya merkezi-fark (HS gibi başlangıçta kabul). φ>0'da
   bölge-tutarlı analitik teğet ŞART (FD bölge sınırında Newton'u kırar — MC dersi).
6. **Doğrula (basitten karmaşığa):** tek-eleman gerilme yolu (triaxial/oedometer) kapalı-formla →
   limit yük (Prandtl/şev) → kuplajlı → sektör-programı benchmark'ı (aşağı §5). <%5 hedefi.

**Tasarım kararı:** etiketli-union (enum+switch) bilinçli — vtable/indirection yok, hot-loop hızlı,
derleyici inline eder. Maliyet: switch iki yerde (plane/axisym). Model sayısı çok artarsa C++20
`concept`-tabanlı kayıt (registry) düşünülebilir; şimdilik etiketli-union sektör kodlarının da yaklaşımı.

## 5. Hedef modeller (bu mimariye oturur)
- **Mohr-Coulomb** ✅ (asal-uzay 4-bölge return, analitik teğet).
- **Hardening Soil** ✅ (büyük oranda; cap+shear çift-yüzey, substepping).
- **HSsmall** ⏳ (HS + küçük-şekildeğiştirme rijitliği G0/γ0.7, Hardin-Drnevich; HS'in üstüne G(γ)).
- **Soft Soil / Modified Cam-Clay** ⏳ (kritik-durum; cap = elips).
- **Tresca / undrained** ✅ (MC, φ=0).
- İleride: Sekiguchi-Ohta, UBCSAND, vb. — hepsi aynı `integrate_point` arayüzünden girer.

İlgili: [[constitutive-models]], [[mohr-coulomb-formulation]], [[hardening-soil-formulation]],
[[multi-code-validation-plan]].
