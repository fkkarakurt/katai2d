"""A pile row's skin springs are tied to the soil around each of them.

The spring between a pile row and the ground is a multiple of the shear modulus of the soil it
lies in. A pile through two layers therefore sheds its load at two rates: where the pile and the
ground both move continuously across the boundary, the skin traction -- the slope of the pile's
axial force -- jumps there by the ratio of the two moduli.

An elastic pile is pushed down through a soft layer into one ten times as stiff, and the slope of
its axial force is read on both sides of the boundary. With one modulus for the whole pile (that
of the soil at its middle, which is what this program did) the two slopes are the same.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def slopes(E_low):
    prj = katai.Project("pile row", element="tri6", mesh_size=0.5, auto_refine=False)
    up = prj.materials.linear_elastic("Upper", E=1.0e4, nu=0.3, gamma=0.0, gamma_sat=0.0)
    low = prj.materials.linear_elastic("Lower", E=E_low, nu=0.3, gamma=0.0, gamma_sat=0.0)
    prj.geometry.polygon([(0, -5), (8, -5), (8, 0), (0, 0)], material=up)
    prj.geometry.polygon([(0, -12), (8, -12), (8, -5), (0, -5)], material=low)
    pile = prj.structures.pile((4, 0), (4, -9), E=3.0e7, diameter=0.4, spacing=2.0, gamma=0.0,
                               name="Pile")
    load = prj.loads.point_load((4, 0), qy=-100.0, name="P")
    prj.initial(procedure="k0", exclude=[pile, load])
    prj.phases.plastic("pile", activate=[pile])
    prj.phases.plastic("load", activate=[load])
    job = katai.Job(prj.build())
    if not job.run():
        return None, job.message()
    sf = [s for s in job.results()[-1].struct_forces if s.name == "Pile"][0]
    st = sorted(((p.y, p.N) for p in sf.stations), reverse=True)

    def slope(y_hi, y_lo):   # least-squares dN/dy over the stations in (y_lo, y_hi)
        pts = [(y, n) for y, n in st if y_lo < y < y_hi]
        my = sum(y for y, _ in pts) / len(pts)
        mn = sum(n for _, n in pts) / len(pts)
        return sum((y - my) * (n - mn) for y, n in pts) / sum((y - my) ** 2 for y, _ in pts)

    return (slope(-3.0, -4.9), slope(-5.1, -7.0)), ""


same, msg = slopes(1.0e4)
check(same is not None, "the pile in one soil runs" + ("" if same else " -- " + msg[:300]))
two, msg = slopes(1.0e5)
check(two is not None, "the pile through two soils runs" + ("" if two else " -- " + msg[:300]))
if same and two:
    r1 = same[1] / same[0]
    r10 = two[1] / two[0]
    check(0.5 < r1 < 1.3,
          f"in one soil the axial force falls at about one rate across the level: ratio {r1:.2f}")
    check(r10 > 4.0 * r1,
          f"into a soil ten times as stiff it falls several times faster below the boundary: "
          f"ratio {r10:.2f} against {r1:.2f} in one soil")

sys.exit(1 if failures else 0)
