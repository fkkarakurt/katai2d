"""A plate that ends on a symmetry axis cannot rotate there.

A slab on elastic ground under a load at its middle is solved whole, and as the half to the right
of its middle with the cut as a side of the model. The half is the whole only if the slab is held
against rotation on the cut: with the rotation free there, its bending moment at the cut is zero
-- at the very section where the whole slab's is largest.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def slab(x0):
    prj = katai.Project("slab", element="tri6", mesh_size=0.4, auto_refine=False)
    soil = prj.materials.linear_elastic("Ground", E=2.0e4, nu=0.3, gamma=0.0, gamma_sat=0.0)
    prj.geometry.polygon([(x0, 0), (8, 0), (8, 4), (x0, 4)], material=soil)   # standard fixities
    plate = prj.structures.plate((max(x0, -3.0), 4), (3, 4), EA=5.0e6, EI=5.0e4, w=0.0, nu=0.0,
                                 name="Slab")
    q = prj.loads.line_load((max(x0, -1.0), 4), (1, 4), qy=-100.0, name="q")
    prj.initial(procedure="k0", exclude=[plate, q])
    prj.phases.plastic("load", activate=[plate, q])
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


def at_middle(job):
    st = [q for s in job.results()[-1].struct_forces for q in s.stations]
    return min(st, key=lambda q: abs(q.x)).M, max(abs(q.M) for q in st)


ok_w, whole = slab(-8.0)
ok_h, half = slab(0.0)
check(ok_w and ok_h, "the whole slab and its half run"
      + ("" if ok_w and ok_h else " -- " + (whole.message() + half.message())[:300]))
if ok_w and ok_h:
    mw, peak_w = at_middle(whole)
    mh, peak_h = at_middle(half)
    check(abs(mw) > 0.9 * peak_w, f"the whole slab's moment is largest at its middle ({mw:.2f} of {peak_w:.2f} kNm/m)")
    check(abs(mh / mw - 1) < 0.03,
          f"the half has that moment on the cut: {mh:.2f} against {mw:.2f} kNm/m "
          f"({100 * (mh / mw - 1):+.2f}%; with the rotation free there it is 0)")

sys.exit(1 if failures else 0)
