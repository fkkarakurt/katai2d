"""The automatic K0 of a model that carries K0nc starts from that K0nc.

A Hardening Soil or Soft Soil material has K0nc as a parameter: it sets the cap, and the
initial stresses of normally consolidated ground lie on that line. The automatic K0 used to be
Jaky's 1 - sin(phi') whatever was entered, so a soil given K0nc = 0.65 at phi' = 30 degrees
started at 0.5 under a cap built for 0.65. Level ground, no water: the ratio of the horizontal
to the vertical stress at any point IS the K0 that was used.
"""
import math
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def ratio(material):
    prj = katai.Project("k0", element="tri6", mesh_size=1.0, auto_refine=False)
    soil = material(prj)
    prj.geometry.polygon([(0, 0), (10, 0), (10, 10), (0, 10)], material=soil,
                         fix=["full", "horizontal", "free", "horizontal"])
    prj.initial(procedure="k0")
    job = katai.Job(prj.build())
    if not job.run():
        return float("nan"), job.message()
    r = job.results()[0]
    # one row of nodal stress per node, (xx, yy, ...): the node nearest mid-depth, mid-width
    n = len(r.node_x)
    i = min(range(n), key=lambda k: abs(r.node_x[k] - 5) + abs(r.node_y[k] - 5))
    row = r.stress[i]
    return float(row[0]) / float(row[1]), ""


hs_kw = dict(E50=20000, Eoed=20000, Eur=60000, m=0.5, c=1.0, phi=30.0, psi=0.0, gamma=18.0,
             nu_ur=0.2, p_ref=100.0, Rf=0.9)
k_entered, msg = ratio(lambda p: p.materials.hardening_soil("HS", k0nc_auto=False, k0nc=0.65,
                                                           **hs_kw))
check(abs(k_entered - 0.65) < 1e-3,
      f"Hardening Soil with K0nc = 0.65 entered starts at K0 = {k_entered:.4f} "
      f"(was 1 - sin 30 = 0.5)" + (" -- " + msg[:200] if msg else ""))
k_auto, msg = ratio(lambda p: p.materials.hardening_soil("HS", **hs_kw))
check(abs(k_auto - (1 - math.sin(math.radians(30)))) < 1e-3,
      f"and with K0nc left automatic at Jaky's: {k_auto:.4f}" + (" -- " + msg[:200] if msg else ""))
k_ss, msg = ratio(lambda p: p.materials.soft_soil("SS", lam_star=0.1, kap_star=0.02, c=1.0, phi=25.0,
                                                  gamma=16.0, k0nc_auto=False, k0nc=0.7))
check(abs(k_ss - 0.7) < 1e-3,
      f"Soft Soil with K0nc = 0.7 entered starts at K0 = {k_ss:.4f}" + (" -- " + msg[:200] if msg else ""))
k_oc, msg = ratio(lambda p: p.materials.hardening_soil("HS", k0nc_auto=False, k0nc=0.65, oc_mode=1,
                                                       OCR=2.0, **hs_kw))
want = 0.65 * 2.0 - 0.2 / 0.8 * 1.0
check(abs(k_oc - want) < 1e-3,
      f"and the overconsolidation terms go on top of it: OCR = 2 gives {k_oc:.4f} "
      f"(K0nc OCR - nu_ur/(1 - nu_ur) (OCR - 1) = {want:.4f})" + (" -- " + msg[:200] if msg else ""))
k_le, msg = ratio(lambda p: p.materials.linear_elastic("LE", E=20000, nu=0.3, gamma=18.0))
check(abs(k_le - 0.5) < 1e-3,
      f"a model with no friction angle takes 0.5: {k_le:.4f}" + (" -- " + msg[:200] if msg else ""))

sys.exit(1 if failures else 0)
