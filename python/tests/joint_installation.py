"""What a joint starts with: gravity loading with a wall in place, and a wall put into wet ground.

  * A model started by GRAVITY LOADING is loaded by gravity whether or not a wall with
    interfaces stands in it. It used to be run by the K0 procedure the moment a joint was active
    in the initial phase, and said nothing: no settlement, and a lateral stress of K0 sigma_v
    where elastic ground under its own weight gives nu / (1 - nu) sigma_v.
  * A wall installed in ground that carries an EXCESS PORE PRESSURE takes over the total stress
    at its line, the effective stress and the pressure. Its joints find the stress their friction
    acts on by taking the pressure off again, so seeded with the effective stress alone they had
    it taken off twice.
  * The joint's own results are readable (interface_forces): they are what both statements are
    checked on.
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


GAMMA, NU, W, H = 18.0, 0.2, 10.0, 6.0
ELASTIC_RATIO = NU / (1.0 - NU)          # 0.25; this material's K0 is 0.5


def ground_ratio(r, at):
    i = min(range(len(r.node_x)), key=lambda k: math.hypot(r.node_x[k] - at[0], r.node_y[k] - at[1]))
    return r.stress[i][0] / r.stress[i][1]


# --- 1. gravity loading with the wall already in place -----------------------------------------
def gravity_start(wall_in_place):
    prj = katai.Project("gravity", element="tri6", mesh_size=1.0, auto_refine=False)
    soil = prj.materials.linear_elastic("Ground", E=2.0e4, nu=NU, gamma=GAMMA, gamma_sat=GAMMA)
    rough = prj.materials.mohr_coulomb("Joint", E=2.0e4, nu=NU, c=3.5, phi=18.0, psi=0.0,
                                       gamma=GAMMA, gamma_sat=GAMMA)
    prj.geometry.polygon([(0, 0), (W, 0), (W, H), (0, H)], material=soil)
    wall = prj.structures.plate((5, H), (5, 2), EA=1.0e7, EI=1.0e5, w=0.0, nu=0.0,
                                interfaces="both", name="Wall", iface_material=rough)
    prj.initial(procedure="gravity", exclude=[] if wall_in_place else [wall])
    job = katai.Job(prj.build())
    job.run()
    return job.results()[0]


free, walled = gravity_start(False), gravity_start(True)
print(f"      gravity loading: settlement {free.max_disp * 1e3:.3f} mm without the wall, "
      f"{walled.max_disp * 1e3:.3f} mm with it; sigma_h / sigma_v far from the wall "
      f"{ground_ratio(free, (0.5, 1.0)):.4f} and {ground_ratio(walled, (0.5, 1.0)):.4f} "
      f"(nu / (1 - nu) = {ELASTIC_RATIO:.4f}, K0 = 0.5)")
check(abs(ground_ratio(free, (0.5, 1.0)) / ELASTIC_RATIO - 1) < 1e-6,
      "without a wall the ground under its own weight has sigma_h = nu / (1 - nu) sigma_v")
check(walled.max_disp > 0.9 * free.max_disp,
      "with the wall in place the ground settles under its weight as well (it did not move at all)")
check(abs(ground_ratio(walled, (0.5, 1.0)) / ELASTIC_RATIO - 1) < 0.10,
      "and far from the wall its lateral stress is the elastic one, not K0 sigma_v")
check(len(walled.interface_forces) > 0 and len(walled.interface_forces[0].stations) > 0,
      "the wall's joints report their stresses (interface_forces)")
if walled.interface_forces:
    top = [q for ir in walled.interface_forces for q in ir.stations if q.y > H - 1e-9]
    check(all(abs(q.sigma_n) < 1e-6 for q in top),
          "a joint that started stress-free carries no normal stress at the ground surface")

# --- 2. a wall installed in ground that carries an excess pore pressure -------------------------
Q = 50.0


def wet_install(drainage):
    prj = katai.Project("wet install", element="tri6", mesh_size=1.0, auto_refine=False)
    soil = prj.materials.linear_elastic("Ground", E=2.0e4, nu=0.3, gamma=GAMMA, gamma_sat=GAMMA,
                                        drainage=drainage)
    rough = prj.materials.mohr_coulomb("Joint", E=2.0e4, nu=0.3, c=3.5, phi=18.0, psi=0.0,
                                       gamma=GAMMA, gamma_sat=GAMMA)
    prj.geometry.polygon([(0, 0), (W, 0), (W, H), (0, H)], material=soil)
    wall = prj.structures.plate((5, H), (5, 2), EA=1.0e7, EI=1.0e5, w=0.0, nu=0.0,
                                interfaces="both", name="Wall", iface_material=rough)
    q = prj.loads.line_load((0, H), (W, H), qy=-Q, name="q")
    prj.initial(procedure="k0", exclude=[wall, q])
    prj.phases.plastic("surcharge", activate=[q])
    prj.phases.plastic("install the wall", activate=[wall])
    job = katai.Job(prj.build())
    job.run()
    loaded, installed = job.results()[1], job.results()[2]
    water = -min(loaded.excess_pore) if len(loaded.excess_pore) else 0.0
    # mid-height of the wall: the joint's transmitted stress against the ground's effective one
    i = min(range(len(loaded.node_x)),
            key=lambda k: math.hypot(loaded.node_x[k] - 5.0, loaded.node_y[k] - 4.0))
    eff_h = -loaded.stress[i][0]
    joint = [-s.sigma_n for ir in installed.interface_forces for s in ir.stations
             if abs(s.y - 4.0) < 1e-6]
    return water, eff_h, joint, installed.max_disp


for drainage in ("drained", "undrained_a"):
    water, eff_h, joint, moved = wet_install(drainage)
    print(f"      {drainage:12s}: water {water:7.3f} kPa, ground sigma'_h {eff_h:7.3f} kPa, joint "
          f"transmits {min(joint):7.3f} .. {max(joint):7.3f} kPa; installing moved {moved * 1e3:.4f} mm")
    # (the wall's head stands under the surcharge and takes its share of it, so installing it is
    # not a nil step here: 0.43 mm, drained and undrained alike)
    check(all(abs(j - (eff_h + water)) < 0.02 * (eff_h + water) for j in joint),
          f"{drainage}: the joint takes over sigma'_h + u = {eff_h + water:.2f} kPa")

sys.exit(1 if failures else 0)
