"""Line contraction of a plate: a shrinkage strain of half the area loss, realised as far as the
surroundings let it.

A plate lies along a weightless strip of ground a million times softer than itself, so the
ground restrains nothing. Two closed statements bound what a contraction C can do:

  * one end held, the other free: the plate shrinks freely -- the free end moves by C/2 x L and
    the plate carries no force;
  * both ends held: nothing can shrink -- the plate does not move and carries N = EA x C/2.

Taking the contraction back in a later phase undoes both.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


L, EA, C = 10.0, 1.0e6, 0.01


def strip(right_held):
    prj = katai.Project("contraction", element="tri6", mesh_size=0.5, auto_refine=False)
    soil = prj.materials.linear_elastic("Soft", E=1.0, nu=0.0, gamma=0.0, gamma_sat=0.0)
    prj.geometry.polygon([(0, 0), (L, 0), (L, 1), (0, 1)], material=soil,
                         fix=["free", "full" if right_held else "free", "free", "full"])
    plate = prj.structures.plate((0, 1), (L, 1), EA=EA, EI=1.0e3, w=0.0, nu=0.0, name="Bar")
    prj.initial(procedure="k0", exclude=[plate])
    prj.phases.plastic("install", activate=[plate])
    prj.phases.plastic("contract", contraction={plate: C})
    prj.phases.plastic("hold")
    prj.phases.plastic("release", contraction={plate: 0.0})
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


def end_ux(res):
    return max((res.displacement[2 * i] for i in range(len(res.node_x))
                if abs(res.node_x[i] - L) < 1e-9 and abs(res.node_y[i] - 1.0) < 1e-9), key=abs)


def normal(res):
    return [q.N for s in res.struct_forces for q in s.stations]


ok_f, free = strip(False)
check(ok_f, "the strip with a free end runs" + ("" if ok_f else " -- " + free.message()[:300]))
if ok_f:
    r = free.results()
    u = end_ux(r[2])
    check(abs(u / (-0.5 * C * L) - 1) < 1e-3,
          f"free to shrink, its end moves by half the contraction times its length: "
          f"{1000 * u:.4f} mm against {-1000 * 0.5 * C * L:.4f}")
    n = max(abs(x) for x in normal(r[2]))
    check(n < 1e-3 * EA * 0.5 * C, f"and it carries no force ({n:.3e} kN/m against EA C/2 = {EA * 0.5 * C:.0f})")
    check(abs(end_ux(r[3])) < 1e-6 * abs(u), "a phase that keeps the contraction moves nothing")
    back = end_ux(r[4])
    check(abs(back / (0.5 * C * L) - 1) < 1e-3,
          f"taking it back moves the end back: {1000 * back:.4f} mm")

ok_h, held = strip(True)
check(ok_h, "the strip held at both ends runs" + ("" if ok_h else " -- " + held.message()[:300]))
if ok_h:
    r = held.results()
    ns = normal(r[2])
    want = EA * 0.5 * C
    check(max(abs(x / want - 1) for x in ns) < 1e-3,
          f"held at both ends it cannot shrink and is in tension by EA C/2: "
          f"{min(ns):.2f} .. {max(ns):.2f} kN/m against {want:.2f}")
    check(max(abs(x) for x in normal(r[4])) < 1e-6 * want,
          "and taking the contraction back takes the force away")

sys.exit(1 if failures else 0)
