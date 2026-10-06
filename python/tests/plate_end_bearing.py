"""End bearing of a plate: no soil plasticity around its lower end.

A plate has no thickness, so a wall pushed down its own axis meets the ground under it only at a
node. With end bearing on, the ground within the wall's equivalent thickness d = sqrt(12 EI/EA) of
its foot answers elastically. Three statements that do not need a tip-resistance formula -- the
option is explicitly not one:

  * in ground that is elastic anyway the option changes nothing at all;
  * in ground that yields, the wall settles less with the option than without;
  * the phase says which zone it made (K2D-A029), and a phase without the wall has none;
  * a wall too thin for the mesh -- no element within its d of the foot -- has no zone: the
    answer is that of the wall without the option, and that is said (K2D-A030).
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def wall(punch, yielding, EI=1.28e6):
    prj = katai.Project("end bearing", element="tri6", mesh_size=0.5, auto_refine=False)
    if yielding:
        soil = prj.materials.mohr_coulomb("Clay", E=1.0e4, nu=0.3, c=8.0, phi=20.0, psi=0.0,
                                          gamma=18.0, gamma_sat=18.0)
    else:
        soil = prj.materials.linear_elastic("Elastic", E=1.0e4, nu=0.3, gamma=18.0, gamma_sat=18.0)
    prj.geometry.polygon([(0, 0), (12, 0), (12, 10), (0, 10)], material=soil)
    # 0.8 m thick: EA = E t, EI = E t^3 / 12  ->  d = 0.8 m
    w = prj.structures.plate((6, 10), (6, 5), EA=2.4e7, EI=EI, w=0.0, nu=0.0, name="Wall",
                             prevent_punching=punch)
    load = prj.loads.point_load((6, 10), qy=-400.0, name="P")
    prj.initial(procedure="k0", exclude=[w, load])
    prj.phases.plastic("wall", activate=[w])
    prj.phases.plastic("load", activate=[load])
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job


def head(job):
    r = job.results()[-1]
    i = min(range(len(r.node_x)), key=lambda k: (r.node_x[k] - 6.0) ** 2 + (r.node_y[k] - 10.0) ** 2)
    return -r.displacement[2 * i + 1]


def codes(res):
    return [d.code for d in res.diagnostics]


ok_e0, e0 = wall(False, False)
ok_e1, e1 = wall(True, False)
check(ok_e0 and ok_e1, "the wall in elastic ground runs with and without end bearing"
      + ("" if ok_e0 and ok_e1 else " -- " + (e0.message() + e1.message())[:300]))
if ok_e0 and ok_e1:
    a, b = head(e0), head(e1)
    check(a > 0 and abs(b / a - 1) < 1e-12,
          f"in elastic ground the option changes nothing: {1000 * a:.6f} and {1000 * b:.6f} mm")

ok_y0, y0 = wall(False, True)
ok_y1, y1 = wall(True, True)
check(ok_y1, "the wall with end bearing carries the load in ground that yields"
      + ("" if ok_y1 else " -- " + y1.message()[:300]))
if ok_y1:
    b = head(y1)
    if ok_y0:
        a = head(y0)
        check(b < 0.97 * a, f"and settles less than without it: {1000 * b:.3f} against {1000 * a:.3f} mm")
    else:
        print(f"      without end bearing the load is not carried at all ({y0.message()[:120]}); "
              f"with it the head settles {1000 * b:.3f} mm")
    r = y1.results()
    check("K2D-A029" in codes(r[1]) and "K2D-A029" in codes(r[2]),
          "the phases with the wall say which zone they made (K2D-A029)")
    check("K2D-A029" not in codes(r[0]), "and the initial phase, without the wall, has none")
if ok_y0:
    check("K2D-A029" not in [c for res in y0.results() for c in codes(res)],
          "a wall without the option makes no zone")

# EA / EI of a 0.1 m section: d = 0.1 m, and no element of a 0.5 m mesh has its centroid that close.
ok_c0, c0 = wall(False, False, EI=2.0e4)
ok_c1, c1 = wall(True, False, EI=2.0e4)
check(ok_c0 and ok_c1, "a wall 0.1 m thick runs with and without the option")
if ok_c0 and ok_c1:
    a, b = head(c0), head(c1)
    cc = [c for res in c1.results() for c in codes(res)]
    check(abs(b / a - 1) < 1e-12 and "K2D-A030" in cc and "K2D-A029" not in cc,
          f"its zone holds no element of this mesh: the answer is the wall's without the option "
          f"({1000 * a:.6f} and {1000 * b:.6f} mm) and that is said (K2D-A030)")

sys.exit(1 if failures else 0)
