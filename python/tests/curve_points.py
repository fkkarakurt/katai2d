"""Points for curves from Python: a point set with ``prj.curve_point`` is in the project, and each
phase's ``result.curve(i)`` gives its displacement at every step.

A weightless linear elastic block under a strip load applied in four steps: the point under the
load moves in proportion to Sum-Mstage (linear elasticity), its curve starts at rest and ends at its
node's displacement in the phase's result, and the curve's step values are the phase's path.
"""
import sys

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


prj = katai.Project("Curves", mesh_size=0.5, auto_refine=False)
soil = prj.materials.linear_elastic("Elastic", E=1.0e4, nu=0.3, gamma=0.0)
prj.geometry.rectangle(0, 0, 10, 5, material=soil)
strip = prj.loads.line_load((0, 5), (2, 5), qy=-100.0, name="Strip")
a = prj.curve_point(1.0, 5.0, name="Under")
b = prj.curve_point(9.0, 2.5)
prj.initial(exclude=[strip])
prj.phases.plastic("Load", activate=[strip], load_steps=4)
pr = prj.build()
check(len(pr.curves) == 2 and pr.curves[0].name == "Under" and pr.curves[1].name == "B",
      "the points are in the project, the second named by its place")
r = prj.run().results()[-1]
curve = r.curve(a)
print("   %d steps; the point under the load ends at uy = %.6e m" % (len(curve) - 1, curve[-1][2] if curve else 0.0))
check(len(curve) == len(r.path_factor) >= 4 and curve[0][1:] == (0.0, 0.0),
      "one entry per step of the path, starting at rest")
end = curve[-1]
node = r.curve_nodes[a]
check(abs(end[1] - r.displacement[2 * node]) <= 1e-12 and abs(end[2] - r.displacement[2 * node + 1]) <= 1e-12,
      "the curve ends at its node's displacement in the phase's result")
worst = max(abs(uy - f * end[2]) for f, ux, uy in curve) / abs(end[2])
check(worst < 1e-9 and end[2] < 0.0, "the point moves in proportion to Sum-Mstage (%.1e)" % worst)
check(r.curve(5) == [], "a point that is not there has no curve")
sys.exit(1 if failures else 0)
