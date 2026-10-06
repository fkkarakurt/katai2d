"""A region can change its material between phases.

Ground that is improved, a fill that sets, a clay that is undrained while it is loaded and
drained in the long term: the phase says which material a region is made of in THAT phase
(``materials={region: material}``), and the region keeps the stresses it carries.

1. A weightless elastic column, laterally confined, loaded twice. The first load meets the soft
   material, the second the stiff one it was changed to: the settlements are q H / E_oed of the
   material that was there when the load came. Changing the material alone moves nothing -- the
   stresses are in equilibrium whatever the stiffness.
2. Changing a region to a material with the SAME parameters changes nothing at all.
3. A hardening material put in place of an elastic one starts from the stresses it finds.
"""
import sys

import katai

failures = 0
last_results = []


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


H, Q1, Q2 = 10.0, 40.0, 60.0
E1, E2, NU = 5000.0, 20000.0, 0.25


def eoed(e):
    return e * (1 - NU) / ((1 + NU) * (1 - 2 * NU))


def column(second, swap_phase=True):
    prj = katai.Project("column", element="tri6", mesh_size=1.0, auto_refine=False)
    soft = prj.materials.linear_elastic("Soft", E=E1, nu=NU, gamma=0.0, gamma_sat=0.0)
    other = second(prj)
    col = prj.geometry.polygon([(0, 0), (2, 0), (2, H), (0, H)], material=soft, name="Column",
                               fix=["full", "horizontal", "free", "horizontal"])
    q1 = prj.loads.line_load((0, H), (2, H), qy=-Q1, name="q1")
    q2 = prj.loads.line_load((0, H), (2, H), qy=-Q2, name="q2")
    prj.initial(procedure="k0", exclude=[q1, q2])
    prj.phases.plastic("First load", activate=[q1])
    if swap_phase:
        prj.phases.plastic("Change the material", materials={col: other})
        prj.phases.plastic("Second load", activate=[q2])
    else:
        prj.phases.plastic("Nothing", activate=[])
        prj.phases.plastic("Second load", activate=[q2])
    job = katai.Job(prj.build())
    ok = job.run()
    res = job.results()
    top = [min(r.displacement[2 * i + 1] for i in range(len(r.node_x)) if abs(r.node_y[i] - H) < 1e-9)
           for r in res] if ok else []
    codes = [[d.code for d in r.diagnostics] for r in res] if ok else []
    global last_results
    last_results = res
    return ok, job.message(), top, codes


ok, msg, top, codes = column(lambda p: p.materials.linear_elastic("Stiff", E=E2, nu=NU, gamma=0.0,
                                                                  gamma_sat=0.0))
check(ok, "the column runs with its material changed between two loads" + ("" if ok else " -- " + msg[:300]))
if ok:
    s1, s2 = Q1 * H / eoed(E1), Q2 * H / eoed(E2)
    check(abs(-top[1] / s1 - 1) < 1e-6,
          f"the first load meets the soft material: {-1000 * top[1]:.4f} mm (q H / E_oed = {1000 * s1:.4f})")
    check(abs(top[2]) < 1e-9,
          f"changing the material alone moves nothing: {1000 * top[2]:.2e} mm")
    check("K2D-A028" in codes[2] and "K2D-A028" not in codes[3],
          "and the phase says the region has another material (K2D-A028), once")
    check(abs(-top[3] / s2 - 1) < 1e-6,
          f"the second load meets the stiff one: {-1000 * top[3]:.4f} mm (q H / E_oed = {1000 * s2:.4f})")
    # ... and the total is the sum of what the phases did, the initial phase counting for nothing.
    tot = [min(r.displacement_total[2 * i + 1] for i in range(len(r.node_x))
               if abs(r.node_y[i] - H) < 1e-9) for r in last_results]
    check(abs(tot[0]) < 1e-12 and abs(tot[1] + s1) < 1e-6 * s1 and abs(tot[3] + s1 + s2) < 1e-6 * (s1 + s2),
          f"the total displacement adds them up: {-1000 * tot[1]:.4f} then {-1000 * tot[3]:.4f} mm "
          f"({1000 * (s1 + s2):.4f})")

ok2, msg2, top2, codes2 = column(lambda p: p.materials.linear_elastic("Soft again", E=E1, nu=NU,
                                                                      gamma=0.0, gamma_sat=0.0))
ok3, msg3, top3, _ = column(lambda p: p.materials.linear_elastic("Unused", E=E2, nu=NU, gamma=0.0,
                                                                 gamma_sat=0.0), swap_phase=False)
check(ok2 and ok3 and abs(top2[3] / top3[3] - 1) < 1e-9,
      "a material with the same parameters in its place changes nothing: "
      + (f"{-1000 * top2[3]:.6f} and {-1000 * top3[3]:.6f} mm" if ok2 and ok3 else (msg2 or msg3)[:300]))

ok4, msg4, top4, _ = column(lambda p: p.materials.hardening_soil(
    "HS", E50=20000, Eoed=20000, Eur=60000, m=0.5, c=5.0, phi=30.0, psi=0.0, gamma=0.0,
    gamma_sat=0.0, nu_ur=0.2, p_ref=100.0, Rf=0.9))
check(ok4 and abs(top4[2]) < 1e-6 and top4[3] < 0.0,
      "a hardening material put in its place starts from the stresses it finds and carries the "
      "second load" + (f" ({1000 * top4[2]:.2e} mm on the change, {-1000 * top4[3]:.3f} mm under the load)"
                       if ok4 else " -- " + msg4[:300]))


# 4. An overconsolidated hardening soil changed to a material of the SAME model keeps its stress
#    history. Re-read from the current stress, the cap would close onto the present state and the
#    overconsolidation would be gone with the change of name: the soil would yield at once under
#    a load its history lets it carry elastically.
def oc_column(swap):
    prj = katai.Project("oc column", element="tri6", mesh_size=1.0, auto_refine=False)
    kw = dict(E50=10000, Eoed=10000, Eur=30000, m=0.5, c=2.0, phi=28.0, psi=0.0, gamma=18.0,
              gamma_sat=18.0, nu_ur=0.2, p_ref=100.0, Rf=0.9, oc_mode=1, OCR=3.0)
    a = prj.materials.hardening_soil("Clay", **kw)
    b = prj.materials.hardening_soil("Clay, renamed", **kw)
    col = prj.geometry.polygon([(0, 0), (2, 0), (2, H), (0, H)], material=a, name="Column",
                               fix=["full", "horizontal", "free", "horizontal"])
    q = prj.loads.line_load((0, H), (2, H), qy=-30.0, name="q")
    prj.initial(procedure="k0", exclude=[q])
    prj.phases.plastic("Load", activate=[q], **({"materials": {col: b}} if swap else {}))
    job = katai.Job(prj.build())
    ok = job.run()
    r = job.results()[-1]
    return ok, job.message(), (min(r.displacement[2 * i + 1] for i in range(len(r.node_x))
                                   if abs(r.node_y[i] - H) < 1e-9) if ok else 0.0)


oka, msga, ua = oc_column(False)
okb, msgb, ub = oc_column(True)
check(oka and okb and abs(ub / ua - 1) < 1e-9,
      "an overconsolidated hardening soil keeps its stress history when it is changed to a "
      "material of the same model: " + (f"{-1000 * ua:.6f} and {-1000 * ub:.6f} mm" if oka and okb
                                        else (msga or msgb)[:300]))

sys.exit(1 if failures else 0)
