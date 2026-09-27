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


def rigid_footing(steps=None):
    prj = katai.Project("footing", kind="axisymmetric", element="tri15", mesh_size=0.4,
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
    return ok, job.message(), force


ok, msg, thrust, pa = active_wall()
check(ok, "c = 0 sand behind a wall moved away carries the whole movement at the default "
          "step count (was: 'collapse load' at 1%)" + ("" if ok else " -- " + msg[:200]))
check(ok and abs(thrust / pa - 1) < 0.05,
      f"and the thrust is Rankine's: {thrust:.2f} against Pa = {pa:.2f} kN/m "
      f"({100 * (thrust / pa - 1):+.2f}%; tri6, 0.5 m)")

ok20, msg20, f20 = rigid_footing()
check(ok20, "a rigid footing pushed 50 mm into sand (15-noded, psi = 0) reaches full settlement "
            "at the default step count (was: out of iterations at 20%)"
      + ("" if ok20 else " -- " + msg20[:200]))
ok60, msg60, f60 = rigid_footing(60)
check(ok60, "and at three times the steps" + ("" if ok60 else " -- " + msg60[:200]))
check(ok20 and ok60 and abs(f20 / f60 - 1) < 0.005,
      f"with the same footing force: {f20:.2f} and {f60:.2f} kN/rad "
      f"({100 * (f20 / f60 - 1):+.3f}%)")



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


# 3. An undrained soil in a Safety run. The search re-solves the ground from the unstressed
#    state, so every trial loaded the undrained clay with its whole self-weight, put that weight
#    into the pore water and reported a factor of safety far too low (0.73 against 1.4 on a
#    published embankment), with no warning. Refused now (K2D-G017); the remedies still run.
ok, codes, _ = slope_safety("undrained_a")
check(not ok and "K2D-G017" in codes,
      "an Undrained (A) soil with phi' > 0 in a Safety phase is refused (K2D-G017)")
ok_d, codes_d, fos_d = slope_safety("drained")
check(ok_d and "K2D-G017" not in codes_d and fos_d > 0,
      f"the drained remedy runs (FoS {fos_d:.3f})")
ok_i, codes_i, fos_i = slope_safety("undrained_a", ignore=True)
check(ok_i and "K2D-G017" not in codes_i and abs(fos_i - fos_d) < 1e-9,
      f"so does ignoring undrained behaviour, to the drained factor bit for bit ({fos_i:.3f})")

sys.exit(1 if failures else 0)
