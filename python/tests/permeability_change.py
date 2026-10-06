"""Change of permeability with void ratio in a consolidation phase: log10(k/k0) = de/ck.

A weightless oedometer column is first compressed drained by Q1, which leaves a uniform
volumetric strain -Q1/E_oed and so a uniform change in void ratio de = (1 + e0)(exp(eps_v) - 1).
A small extra load Q2 is then consolidated. With a change of permeability the column has, in
that phase, the permeability k0 * 10^(de/ck) everywhere -- so its settlement after a time t must
be that of the same column WITHOUT a change of permeability after the time t * 10^(de/ck). The two
runs take the same number of time steps, so they solve the same dimensionless problem.
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


H, E, Q1, Q2, E0, CK, K = 6.0, 10000.0, 200.0, 1.0, 0.5, 0.1, 0.002
DE = (1.0 + E0) * math.expm1(-Q1 / E)
KFAC = 10.0 ** (DE / CK)
T = 0.5 * H * H / (K * KFAC * E / 9.81)      # Tv = 0.5 with the reduced permeability


def column(ck, duration, model="mc"):
    prj = katai.Project("ck", element="tri6", mesh_size=0.4, auto_refine=False)
    kw = dict(E=E, nu=0.0, gamma=0.0, gamma_sat=0.0, k=K, e_init=E0)
    if ck:
        kw["ck"] = ck
    if model == "mc":   # strong enough to stay elastic: the path with a state, an elastic answer
        soil = prj.materials.mohr_coulomb("Soil", c=1.0e5, phi=30.0, psi=0.0, **kw)
    else:
        soil = prj.materials.linear_elastic("Soil", **kw)
    prj.geometry.rectangle(0.0, 0.0, 1.0, H, material=soil, name="Column",
                           flow={"top": ("head", H), "bottom": "closed",
                                 "right": "closed", "left": "closed"})
    q1 = prj.loads.line_load((0.0, H), (1.0, H), qy=-Q1, name="Preload")
    q2 = prj.loads.line_load((0.0, H), (1.0, H), qy=-Q2, name="Extra")
    prj.initial(procedure="k0", exclude=[q1, q2])
    prj.phases.plastic("Preload", activate=[q1])
    prj.phases.consolidation("Consolidate", duration=duration, steps=40, activate=[q2],
                             load_at_start=True)
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


ok_a, a = column(CK, T)
ok_b, b = column(0.0, T * KFAC)
ok_c, c = column(0.0, T)
check(ok_a and ok_b and ok_c, "the three columns run"
      + ("" if ok_a and ok_b and ok_c else " -- " + (a.message() + b.message() + c.message())[:300]))
if ok_a and ok_b and ok_c:
    sa, sb, sc = (list(j.results()[-1].consol_settlement) for j in (a, b, c))
    s_inf = Q2 * H / E
    worst = max(abs(x - y) for x, y in zip(sa, sb)) / s_inf
    check(len(sa) == len(sb) and worst < 5e-3,
          f"with a change of permeability the settlement after t is that of the constant column "
          f"after t x {KFAC:.4f}: the two series differ by at most {100 * worst:.3f}% of the final "
          f"settlement (U at the end {sa[-1] / s_inf:.4f} and {sb[-1] / s_inf:.4f})")
    mid = len(sa) // 4
    check(sc[mid] / sa[mid] > 1.15,
          f"and it is slower than the constant column at the same time: a quarter of the way "
          f"U = {sa[mid] / s_inf:.4f} against {sc[mid] / s_inf:.4f}")
    codes = [d.code for r in a.results() for d in r.diagnostics]
    check("K2D-M007" not in codes, "a soil whose void ratio is followed is not warned about")

ok_d, d = column(CK, T, model="le")
codes_d = [x.code for r in d.results() for x in r.diagnostics] if ok_d else []
check(ok_d and "K2D-M007" in codes_d,
      "entered where every soil is linear elastic it is reported as not applied (K2D-M007)"
      + ("" if ok_d else " -- " + d.message()[:300]))

sys.exit(1 if failures else 0)
