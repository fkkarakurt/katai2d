"""Three answers that were wrong while nothing in the model was, from the script.

1. A cohesionless soil reported as collapsing when nothing had failed. With c = 0 the
   Mohr-Coulomb apex is the stress-free state and its tangent is exactly zero, so a node
   whose stress points all lost their stress had no stiffness; the linear solver refused,
   three halvings reached the floor, and the phase printed "the incremental limit
   (collapse) load" at 1% of a 4 mm wall movement. The case is a smooth rigid wall moved
   away from c = 0, phi = 41 sand, for which Rankine gives the answer in closed form:
   Pa = 1/2 Ka gamma H^2, Ka = (1 - sin phi)/(1 + sin phi).

2. A rigid footing pushed into sand that ran out of iterations. A prescribed displacement
   entered each increment through the fixed nodes alone, so the first iterate strained
   one row of elements by the whole increment -- an out-of-balance force 22 times the
   reference -- and Newton spent its budget walking out of stresses the converged field
   never has. At 20 load steps the phase stopped at 20% of the settlement, at 100 steps at
   37%. There is no closed form for the force at 50 mm; what is asserted is the property
   an answer has and a patience setting does not: it does not move with the step count.

3. An undrained soil in a Safety run, reported far too weak (see below; K2D-G017).
"""
import math
import sys

import katai

failures = 0
last_codes = []


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def active_wall(steps=None):
    phi, gamma, H, W = 41.0, 20.4, 8.3, 20.0
    prj = katai.Project("active", element="tri6", mesh_size=0.5, auto_refine=False)
    sand = prj.materials.mohr_coulomb("Sand", E=5e4, nu=0.3, c=0.0, phi=phi, psi=phi,
                                      gamma=gamma, k0_auto=False, k0=0.5)
    prj.geometry.polygon([(0, 0), (W, 0), (W, H), (0, H)], material=sand,
                         fix=["full", "horizontal", "free", "free"])
    move = prj.displacements.line((0, 0), (0, H), ux=-0.004, name="Wall")
    rest = prj.displacements.line((0, 0), (0, H), ux=0.0, name="Wall at rest")
    prj.initial(procedure="k0", exclude=[move, rest])
    prj.phases.plastic("At rest", activate=[rest])
    prj.phases.plastic("Move away", activate=[move], deactivate=[rest],
                       **({"load_steps": steps} if steps else {}))
    job = katai.Job(prj.build())
    ok = job.run()
    r = job.results()[-1]
    thrust = abs(sum(r.reaction[2 * i] for i in range(len(r.node_x)) if abs(r.node_x[i]) < 1e-9)) \
        if ok else float("nan")
    ka = (1 - math.sin(math.radians(phi))) / (1 + math.sin(math.radians(phi)))
    return ok, job.message(), thrust, 0.5 * ka * gamma * H * H


def rigid_footing(steps=None, mesh=0.4, element="tri15"):
    prj = katai.Project("footing", kind="axisymmetric", element=element, mesh_size=mesh,
                        auto_refine=True)
    sand = prj.materials.mohr_coulomb("Sand", E=13000, nu=0.3, c=1.0, phi=30.0, psi=0.0,
                                      gamma=17.0, gamma_sat=20.0)
    prj.geometry.polygon([(0, 0), (5, 0), (5, 4), (0, 4)], material=sand,
                         fix=["full", "horizontal", "free", "horizontal"])
    prj.water.table(2.0)
    push = prj.displacements.line((0, 4), (1, 4), ux=0.0, uy=-0.05, name="Footing")
    prj.initial(procedure="k0", exclude=[push])
    prj.phases.plastic("Push", activate=[push], **({"load_steps": steps} if steps else {}))
    job = katai.Job(prj.build())
    ok = job.run()
    r = job.results()[-1]
    force = -sum(r.reaction[2 * i + 1] for i in range(len(r.node_x))
                 if abs(r.node_y[i] - 4) < 1e-9 and r.node_x[i] <= 1 + 1e-9) if ok else float("nan")
    global last_codes
    last_codes = [d.code for res in job.results() for d in res.diagnostics]
    return ok, job.message(), force


ok, msg, thrust, pa = active_wall()
check(ok, "c = 0 sand behind a wall moved away carries the whole movement at the default "
          "step count (was: 'collapse load' at 1%)" + ("" if ok else " -- " + msg[:200]))
check(ok and abs(thrust / pa - 1) < 0.05,
      f"and the thrust is Rankine's: {thrust:.2f} against Pa = {pa:.2f} kN/m "
      f"({100 * (thrust / pa - 1):+.2f}%; tri6, 0.5 m)")

ok20, msg20, f20 = rigid_footing()
codes15 = last_codes
check(ok20, "a rigid footing pushed 50 mm into sand (15-noded, psi = 0) reaches full settlement "
            "at the default step count (was: out of iterations at 20%)"
      + ("" if ok20 else " -- " + msg20[:200]))
ok60, msg60, f60 = rigid_footing(60)
check(ok60, "and at three times the steps" + ("" if ok60 else " -- " + msg60[:200]))
check(ok20 and ok60 and abs(f20 / f60 - 1) < 0.005,
      f"with the same footing force: {f20:.2f} and {f60:.2f} kN/rad "
      f"({100 * (f20 / f60 - 1):+.3f}%)")
# The same footing on a finer mesh. An increment Newton could not close was handed to an
# elastic-stiffness iteration before it was made smaller; under an imposed displacement that
# sent the phase down another path and it stopped at 13% of the settlement -- on this mesh only,
# the 0.4 m one above went through. A finer mesh lowers the force, by a few percent.
okf, msgf, ff = rigid_footing(mesh=0.25)
check(okf, "and on a 0.25 m mesh it reaches full settlement too (was: stopped at 13%)"
      + ("" if okf else " -- " + msgf[:200]))
check(ok20 and okf and 0.95 < ff / f20 < 1.0,
      f"with a force a few percent under the coarser mesh's: {ff:.2f} against {f20:.2f} kN/rad "
      f"({100 * (ff / f20 - 1):+.2f}%)")
# The 6-noded element on the same footing: too stiff in axisymmetry once the sand flows at
# constant volume, so its force is higher than the 15-noded one's -- the unsafe side -- and the
# run has to say so (K2D-A027). The 15-noded run must not carry the warning.
ok6, msg6, f6 = rigid_footing(element="tri6")
check(ok6 and "K2D-A027" in last_codes and "K2D-A027" not in codes15,
      "6-noded elements in an axisymmetric model with a yielding soil are warned about "
      "(K2D-A027), 15-noded ones are not" + ("" if ok6 else " -- " + msg6[:200]))
check(ok6 and ok20 and f6 > 1.03 * f20,
      f"because their footing force is the higher one: {f6:.2f} against {f20:.2f} kN/rad "
      f"({100 * (f6 / f20 - 1):+.2f}%)")



def slope_safety(drainage, *, ignore=False):
    """A clay slope with a Safety phase; the clay's drainage type is the variable."""
    prj = katai.Project("slope", element="tri6", mesh_size=1.0, auto_refine=False)
    kw = dict(E=5600, nu=0.3, c=10.0, phi=25.0, gamma=13.0, gamma_sat=13.0, drainage=drainage)
    clay = prj.materials.mohr_coulomb("Clay", **kw)
    prj.geometry.polygon([(0, 0), (30, 0), (30, 6), (18, 6), (10, 2), (0, 2)], material=clay,
                         fix=["full", "horizontal", "free", "free", "free", "horizontal"])
    prj.initial(procedure="k0")
    prj.phases.safety("Safety", **({"ignore_undrained": True} if ignore else {}))
    job = katai.Job(prj.build())
    ok = job.run()
    codes = [d.code for r in job.results() for d in r.diagnostics]
    return ok, codes, job.results()[-1].fos


# 3. An undrained soil in a Safety run. The search used to re-solve the ground from the
#    unstressed state, so every trial loaded the undrained clay with its whole self-weight, put
#    that weight into the pore water and reported a factor of safety far too low (0.73 against
#    1.4 on a published embankment), with no warning; it was refused (K2D-G017) until the search
#    could start from the phase before. A Safety phase now starts from that phase's state, so the
#    undrained soil runs and generates only what the strength reduction itself makes it generate.
ok_u, codes_u, fos_u = slope_safety("undrained_a")
check(ok_u and "K2D-G017" not in codes_u and fos_u > 0,
      f"an Undrained (A) soil with phi' > 0 in a Safety phase after the initial phase runs "
      f"from that phase (FoS {fos_u:.3f})")
ok_d, codes_d, fos_d = slope_safety("drained")
check(ok_d and "K2D-G017" not in codes_d and fos_d > 0,
      f"the drained slope runs (FoS {fos_d:.3f})")
ok_i, codes_i, fos_i = slope_safety("undrained_a", ignore=True)
check(ok_i and "K2D-G017" not in codes_i and abs(fos_i - fos_d) < 1e-9,
      f"and ignoring undrained behaviour gives the drained factor bit for bit ({fos_i:.3f}) -- "
      f"the initial phase is drained whatever the drainage type, so the two start alike")

sys.exit(1 if failures else 0)
