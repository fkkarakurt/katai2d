"""Consolidation of a soil with small-strain stiffness, through a reversal.

A strip load is placed on a sand with small-strain stiffness and then half of it is taken off
again. The unloading is run twice: as a drained Plastic phase, and as a Consolidation phase long
enough for every excess pore pressure to go. The second is the first with time in it, so where it
ends it has to have moved the ground as the first did.

It is the reversal that makes this a test: the stiffness of a step follows the direction of the
step's strain increment against the strain history, and a time step that stalls on a stiffness
following its own iterate is continued with the stiffness read from a fixed increment -- the same
device as in a Plastic phase.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def run(consolidate):
    prj = katai.Project("small strain consolidation", element="tri6", mesh_size=1.0, auto_refine=False)
    sand = prj.materials.hs_small("Sand", E50=25e3, Eoed=25e3, Eur=75e3, G0=80e3, gamma07=0.2e-3,
                                  c=1.0, phi=31, psi=1, m=0.5, nu_ur=0.2, p_ref=100, Rf=0.9,
                                  gamma=17, gamma_sat=20, k=1.0, e_init=0.6)
    prj.geometry.rectangle(0.0, 0.0, 8.0, 6.0, material=sand, name="Ground",
                           flow={"top": ("head", 6.0), "bottom": "closed", "right": "closed",
                                 "left": "closed"})
    prj.water.table(6.0)
    q1 = prj.loads.line_load((0, 6), (2, 6), qy=-40.0, name="q1")
    q2 = prj.loads.line_load((0, 6), (2, 6), qy=+20.0, name="q2")   # half of it taken off again
    prj.initial(procedure="k0", exclude=[q1, q2])
    prj.phases.plastic("load", activate=[q1])
    if consolidate:
        prj.phases.consolidation("unload", duration=50.0, steps=20, activate=[q2], load_at_start=True)
    else:
        prj.phases.plastic("unload", activate=[q2])
    job = katai.Job(prj.build())
    return job.run(), job


def heave(job):
    r = job.results()[-1]
    i = min(range(len(r.node_x)), key=lambda k: r.node_x[k] ** 2 + (r.node_y[k] - 6.0) ** 2)
    return r.displacement[2 * i + 1]


ok_p, p = run(False)
ok_c, c = run(True)
check(ok_p, "the drained unloading runs" + ("" if ok_p else " -- " + p.message()[:300]))
check(ok_c, "the consolidated unloading runs" + ("" if ok_c else " -- " + c.message()[:300]))
if ok_p and ok_c:
    hp, hc = heave(p), heave(c)
    check(hp > 0.0 and abs(hc / hp - 1.0) < 0.03,
          f"consolidated to rest the ground has heaved as in the drained phase: "
          f"{1000 * hc:.3f} mm against {1000 * hp:.3f} mm")
    rc = c.results()[-1]
    left = max(abs(x) for x in rc.excess_pore) if len(rc.excess_pore) else 0.0
    check(left < 0.05, f"and no excess pore pressure is left: {left:.4f} kPa")

sys.exit(1 if failures else 0)
