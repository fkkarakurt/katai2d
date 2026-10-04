"""A phase inherits what the phase before it switched off.

The script surface promises that activation is inherited: a phase lists only what it
changes. It wrote a phase's activation vector only when that phase toggled something of
that class -- and to the engine an EMPTY vector means "everything active", not "as
before". So the first phase after an excavation that did not mention regions brought the
excavated soil back: on this linear-elastic block the floor heaved 13.3 mm when the top
3 m came out, and the next phase, which changed nothing, settled it 12.9 mm back. An
anchor not yet installed was switched on by any phase that did not mention structures,
prestress and all.
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


prj = katai.Project("inheritance", element="tri6", mesh_size=1.0, auto_refine=False)
soil = prj.materials.linear_elastic("Soil", E=5e4, nu=0.3, gamma=19.6, gamma_sat=19.6,
                                    k0_auto=False, k0=0.5)
F, H, X = "free", "horizontal", "full"
top = prj.geometry.polygon([(0, -3), (10, -3), (10, 0), (0, 0)], material=soil, name="Top",
                           fix=[F, F, F, H])
prj.geometry.polygon([(10, -3), (30, -3), (30, 0), (10, 0)], material=soil, name="Top right",
                     fix=[F, H, F, F])
prj.geometry.polygon([(0, -20), (30, -20), (30, -3), (10, -3), (0, -3)], material=soil,
                     name="Base", fix=[X, H, F, F, H])
strut = prj.structures.anchor((10, -1), (20, -6), EA=3e5, spacing=2, prestress=100, name="Tie")
prj.initial(procedure="k0", exclude=[strut])
prj.phases.plastic("Excavate", deactivate=[top])
prj.phases.plastic("Nothing changes")
prj.phases.plastic("Install the tie", activate=[strut])

pr = prj.build()
check(list(pr.phases[1].poly_active) == [0, 1, 1],
      "the phase after the excavation still has the excavated region switched off")
check(list(pr.phases[0].struct_active) == [0] and list(pr.phases[1].struct_active) == [0],
      "and the tie stays out until the phase that installs it")

job = katai.Job(pr)
check(job.run(), "the staged run completes")
res = job.results()
x, y = res[1].node_x, res[1].node_y
floor = min(range(len(x)), key=lambda i: math.hypot(x[i] - 5, y[i] + 3))
heave = 1000 * res[1].displacement[2 * floor + 1]
check(heave > 5.0, f"the excavation heaves its floor ({heave:.3f} mm)")
check(res[2].max_disp == 0.0,
      f"a phase that changes nothing moves nothing (max|u| = {1000 * res[2].max_disp:.6f} mm)")

# THE WATER CARRIES OVER TOO. A phase that names no water keeps the water of the phase before it,
# as a new Studio phase does. It used to fall back to the project's table: a Safety phase after a
# dewatering step filled the excavation again and applied that as a change of its own.
w = katai.Project("water", element="tri6", mesh_size=1.0, auto_refine=False)
ws = w.materials.linear_elastic("Soil", E=5e4, nu=0.3, gamma=18.0, gamma_sat=20.0,
                                k0_auto=False, k0=0.5)
w.geometry.polygon([(0, -20), (30, -20), (30, 0), (0, 0)], material=ws, name="Block",
                   fix=[X, H, F, H])
w.water.table(-2.0)
w.initial(procedure="k0")
w.phases.plastic("Dewater", water=-6.0)
w.phases.plastic("Nothing changes")
w.phases.plastic("Back to the project's table", water="project")
wp = w.build()
check(wp.phases[1].water_override and list(wp.phases[1].wy) == [-6.0, -6.0],
      "a phase that names no water keeps the dewatered table")
check(not wp.phases[2].water_override, 'water="project" returns to the project table')
wj = katai.Job(wp)
check(wj.run(), "the dewatering run completes")
wr = wj.results()
check(wr[1].max_disp > 1e-4, f"lowering the table moves the ground ({1000 * wr[1].max_disp:.3f} mm)")
check(wr[2].max_disp == 0.0,
      f"the phase after it, naming no water, moves nothing (max|u| = {1000 * wr[2].max_disp:.6f} mm)")

sys.exit(1 if failures else 0)
