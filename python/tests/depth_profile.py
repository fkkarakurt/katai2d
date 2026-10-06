"""A stiffness that increases with depth is the reference value above the reference level.

The increment is per unit of depth BELOW y_ref; above it the stiffness is E_ref. A weightless,
laterally confined elastic column 10 m tall with its reference level at mid-height has the
closed form: the upper half compresses at E_ref, the lower half at E_ref + E_inc d,
    s = q/f [ 5/E_ref + ln(1 + 5 E_inc/E_ref)/E_inc ],   f = E_oed/E = (1-nu)/((1+nu)(1-2nu)).
Extrapolated upwards -- as it used to be -- the upper half would be softer than E_ref and the
settlement larger.
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


H, Q, EREF, EINC, NU = 10.0, 100.0, 4000.0, 1200.0, 0.25
prj = katai.Project("profile", element="tri15", mesh_size=0.5, auto_refine=False)
soil = prj.materials.linear_elastic("Soil", E=EREF, nu=NU, gamma=0.0, gamma_sat=0.0,
                                    E_inc=EINC, y_ref=5.0)
prj.geometry.polygon([(0, 0), (1, 0), (1, H), (0, H)], material=soil,
                     fix=["full", "horizontal", "free", "horizontal"])
q = prj.loads.line_load((0, H), (1, H), qy=-Q, name="q")
prj.initial(procedure="k0", exclude=[q])
prj.phases.plastic("Load", activate=[q])
job = katai.Job(prj.build())
ok = job.run()
check(ok, "the column with a stiffness increasing below mid-height runs" + ("" if ok else " -- " + job.message()[:300]))
if ok:
    r = job.results()[-1]
    s = -min(r.displacement[2 * i + 1] for i in range(len(r.node_x)) if abs(r.node_y[i] - H) < 1e-9)
    f = (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    want = Q / f * (5.0 / EREF + math.log(1 + 5.0 * EINC / EREF) / EINC)
    old = Q / f * (math.log((EREF + 5.0 * EINC) / max(EREF - 5.0 * EINC, 1e-9)) / EINC
                   if EREF > 5.0 * EINC else float("inf"))
    check(abs(s / want - 1) < 2e-3,
          f"its settlement is that of E_ref above the reference level and E_ref + E_inc d below: "
          f"{1000 * s:.3f} mm against {1000 * want:.3f} ({100 * (s / want - 1):+.3f}%)")
    check("K2D-M006" not in [d.code for res in job.results() for d in res.diagnostics],
          "and a model that reads the increment says nothing about it")

# A hardening soil does not read a stiffness increment: its stiffness is the stress-dependent
# law. The box is not dropped in silence.
prj2 = katai.Project("unread", element="tri6", mesh_size=1.0, auto_refine=False)
hs = prj2.materials.hardening_soil("HS", E50=20000, Eoed=20000, Eur=60000, m=0.5, c=1.0, phi=30.0,
                                   psi=0.0, gamma=18.0, nu_ur=0.2, p_ref=100.0, Rf=0.9, E_inc=500.0,
                                   y_ref=2.0)
prj2.geometry.polygon([(0, 0), (4, 0), (4, 4), (0, 4)], material=hs,
                      fix=["full", "horizontal", "free", "horizontal"])
prj2.initial(procedure="k0")
job2 = katai.Job(prj2.build())
ok2 = job2.run()
codes2 = [d.code for res in job2.results() for d in res.diagnostics] if ok2 else []
check(ok2 and "K2D-M006" in codes2,
      "an increment entered on a model that does not read it is reported (K2D-M006)"
      + ("" if ok2 else " -- " + job2.message()[:300]))

sys.exit(1 if failures else 0)
