"""Soft Soil with Undrained (A): the pore fluid of a soil whose stiffness follows the stress.

The model has no Young's modulus -- its bulk modulus is (p' + c cot phi)/kappa* -- so the pore
fluid cannot be sized from one, and the combination used to be refused. The fluid now follows
the skeleton, stress point by stress point. A laterally confined column under a uniform
surcharge has the answers:

1. loaded undrained, the water takes the load: the excess pore pressure is the surcharge, all
   the way down, and the column barely settles;
2. left to consolidate, it ends where the same column ends when it is loaded drained.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


H, Q = 10.0, 20.0


def column(drainage, consolidate):
    prj = katai.Project("column", element="tri6", mesh_size=0.5, auto_refine=False)
    clay = prj.materials.soft_soil("Clay", lam_star=0.10, kap_star=0.02, c=1.0, phi=25.0,
                                   gamma=16.0, gamma_sat=16.0, k=1e-3, drainage=drainage)
    prj.geometry.polygon([(0, 0), (1, 0), (1, H), (0, H)], material=clay,
                         fix=["full", "horizontal", "free", "horizontal"],
                         flow=["closed", "closed", ("head", H), "closed"])
    q = prj.loads.line_load((0, H), (1, H), qy=-Q, name="q")
    prj.water.table(H)
    prj.initial(procedure="k0", exclude=[q])
    prj.phases.plastic("Load", activate=[q])
    if consolidate:
        prj.phases.consolidation("Consolidate", until_excess_pore=0.2)
    job = katai.Job(prj.build())
    ok = job.run()
    res = job.results() if ok else []
    top = [min(r.displacement[2 * i + 1] for i in range(len(r.node_x)) if abs(r.node_y[i] - H) < 1e-9)
           for r in res]
    return ok, job.message(), res, top


oku, msgu, resu, topu = column("undrained_a", True)
check(oku, "an Undrained (A) Soft Soil column is loaded and consolidated"
      + ("" if oku else " -- " + msgu[:400]))
okd, msgd, resd, topd = column("drained", False)
check(okd, "and the drained column is loaded" + ("" if okd else " -- " + msgd[:400]))
if oku and okd:
    # The nodal stress is the EFFECTIVE stress. If the water takes the surcharge, the vertical
    # effective stress at mid-depth is what it was before the load; drained, it rises by q.
    def syy_mid(r):
        n = len(r.node_x)
        i = min(range(n), key=lambda k: abs(r.node_x[k] - 0.5) + abs(r.node_y[k] - 0.5 * H))
        return float(r.stress[i][1])
    du = syy_mid(resu[1]) - syy_mid(resu[0])
    dd = syy_mid(resd[1]) - syy_mid(resd[0])
    check(abs(du) < 0.03 * Q and abs(dd + Q) < 0.03 * Q,
          f"loaded undrained, the water takes the surcharge: the vertical effective stress at "
          f"mid-depth changes by {du:+.3f} kPa, against {dd:+.3f} drained (q = {Q:.0f})")
    check(abs(topu[1]) < 0.05 * abs(topd[1]),
          f"and the column barely settles: {-1000 * topu[1]:.2f} mm against {-1000 * topd[1]:.2f} "
          f"drained")
    total = topu[1] + topu[2]
    check(abs(total / topd[1] - 1) < 0.03,
          f"consolidated, it ends where the drained column does: {-1000 * total:.2f} against "
          f"{-1000 * topd[1]:.2f} mm ({100 * (total / topd[1] - 1):+.2f}%)")

sys.exit(1 if failures else 0)
