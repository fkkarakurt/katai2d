# Embedded Beam (Gömülü Kiriş / Pile Row) — Formülasyon (Faz A.4)

Kazık / kaya bulonu / zemin çivisi: bir kiriş elemanının çevre zeminle **keyfi yönelimde, mesh-uyumsuz**
etkileşimi. PLAXIS'in "embedded beam row" elemanı. Etkileşim **skin** (gövde sürtünmesi) ve **foot** (uç
taşıma) gömülü-interface elemanlarıyla. İnsanların kontrol ettiği özellik (kazık grupları).

**Kaynak (kilitli):** PLAXIS 2D 2025.1 Scientific Manual §7.5 (Eq 7-51…7-64; Sadek & Shahrour 2004);
Reference Manual §5.6.3 (geometri) + §6.6 (malzeme: skin/foot bearing capacity, interface stiffness factor);
Sluis (2012). Konvansiyon tension-poz.

## 1. Temel fikir (Sci §7.5)
Kiriş hacim elemanını keyfi konum/yönelimde KESER → eleman içine 3 ek "sanal" düğüm girer. Kiriş Timoshenko
kirişi (mevcut plate çekirdeği). Etkileşim: kiriş düğümleri ↔ zeminin KİRİŞ KONUMUNDAKİ değeri (mesh düğümü
DEĞİL). Bu yüzden zemin şekil fonksiyonları N_s kirişin geçtiği noktada değerlendirilir (non-conforming).

## 2. Kinematik (Eq 7-51…7-56)
u_s = N_s v_s (zemin), u_b = N_b v_b (kiriş). Göreli deplasman:
Δu_rel = Δu_b − Δu_s = N_b Δv_b − N_s Δv_s = N_rel Δv_rel,  N_rel = [N_b  −N_s],  Δv_rel = [Δv_b; Δv_s].

## 3. Skin (gövde) etkileşimi (Eq 7-52…7-58)
Δt_skin = T_skin Δu_rel (T_skin = global skin rijitlik). Sanal iş →
K_skin = ∫ N_rel^T T_skin N_rel dS, dört alt-blok (Eq 7-58):
K_bb=∫N_b^T T N_b, K_bs=∫N_b^T T N_s, K_sb=∫N_s^T T N_b, K_ss=∫N_s^T T N_s.
**Newton-Cotes integrasyon** (Tablo 7-6; bizim interface ile AYNI: 3-düğüm 1/3,4/3,1/3; 5-düğüm
7/45,32/45,12/45) — düğüm-çiftleri ayrık. Plastisitede yalnız elastik kısım rijitlikte; plastisite iteratif.
Skin kapasitesi T_max (Ref §6.6.3) — eksenel sürtünme limiti.

## 4. Foot (uç) etkileşimi (Eq 7-60…7-64)
Δt_foot = D_foot Δu_rel (uç yay rijitliği). K_foot = N_rel^T D_foot N_rel, dört alt-blok (Eq 7-64; işaretler
K_bs=−N_b^T D N_s vb). Foot kapasitesi F_max (uç taşıma).

## 5. KATAI uygulama gereksinimleri (yeni altyapı)
Mevcut plate/interface mesh-UYUMLU (tri6/tri15 kenarı). Embedded beam mesh-UYUMSUZ → GEREKLİ:
1. **Nokta-konumlandırma:** kiriş çizgisi hangi soil elemanından geçer (eleman arama).
2. **Ters haritalama** (fiziksel (x,y) → yerel (ξ,η)) tri6/tri15 için (eğri eleman → Newton iterasyonu);
   sonra N_s(ξ,η) sanal-düğüm interpolasyonu.
3. **Embedded interface eleman:** kiriş düğümü ↔ soil-eleman-düğümleri (N_s ile), skin T_skin + foot D_foot,
   Newton-Cotes, T_max/F_max plastisite (interface return-mapping'e benzer).
4. Kiriş = mevcut plate çekirdeği (3/5-düğüm Timoshenko) ek-DOF'larda.
5. Solver assembly: K_skin/K_foot blokları (kiriş ek-DOF + soil mesh DOF karışık).

## 6. Doğrulama planı
- Eksenel yüklü tek kazık: uç + gövde direnci → yük-oturma; T_max/F_max kapasitesi.
- Yanal yüklü kazık (p-y benzeri) → sehim/moment.
- PLAXIS Tutorial / Ref §6.6 örnek parametreleriyle kıyas (bant/orta).

## 7. Karmaşıklık notu
Non-conforming kuplaj (nokta-konumlandırma + ters-haritalama) yeni geometri altyapısı ister; mesh-uyumlu
plate/interface'ten belirgin daha büyük. Aceleye getirilmemeli (kullanıcı: "hata istemiyorum") → ayrı,
dikkatli, adım-adım doğrulanan bir parça. Önce ters-haritalama helper'ı izole test edilmeli.

İlgili: [[structural-plate-formulation]] (plate çekirdeği), [[interface-formulation]] (Newton-Cotes, return),
[[plaxis-gap-analysis]] Faz A.4, [[material-model-architecture]].
