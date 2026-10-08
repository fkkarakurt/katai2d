"""A support beside undrained ground carries the water's share too.

The stress a support answers is the total one: the effective stress of the skeleton and the excess
pore pressure the ground carries. Three statements, each with an answer that does not come from
this program:

  * a laterally confined column loaded undrained puts the same force on its base as the same
    column loaded drained -- the whole surcharge -- although its effective stress barely moves;
  * a strip footing pushed into undrained clay (s_u, weightless) reaches Prandtl's
    (2 + pi) s_u B, the same force the model reaches when the clay is declared drained with the
    same strength;
  * a joint slips on the EFFECTIVE normal stress. A block held between rigid sides on a rough
    base is loaded and then pushed sideways: loaded drained, the base friction rises by
    q tan(phi_i); loaded undrained, the water takes the load and the friction is that of the
    block with no load on it -- although the stress the joint transmits has risen by all of q.
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


# --- 1. the confined column -------------------------------------------------------------------
GAMMA, Q, W, H = 18.0, 50.0, 4.0, 1.0


def column_base(drainage):
    prj = katai.Project("column", element="tri6", mesh_size=0.5, auto_refine=False)
    soil = prj.materials.linear_elastic("Ground", E=2.0e4, nu=0.3, gamma=GAMMA, gamma_sat=GAMMA,
                                        drainage=drainage)
    prj.geometry.polygon([(0, 0), (W, 0), (W, H), (0, H)], material=soil,
                         fix=["full", "horizontal", "free", "horizontal"])
    q = prj.loads.line_load((0, H), (W, H), qy=-Q, name="q")
    prj.initial(procedure="k0", exclude=[q])
    prj.phases.plastic("surcharge", activate=[q])
    job = katai.Job(prj.build())
    job.run()
    before, after = job.results()[0], job.results()[-1]
    base = [i for i in range(len(after.node_y)) if after.node_y[i] < 1e-9]
    rise = (sum(after.reaction[2 * i + 1] for i in base) -
            sum(before.reaction[2 * i + 1] for i in base))
    # what the water took and what the skeleton took, node by node (both as pressures)
    water = [-p for p in after.excess_pore]
    skeleton = [-(after.stress[i][1] - before.stress[i][1]) for i in range(len(after.node_y))]
    return rise, water, skeleton


(rise_d, water_d, _), (rise_u, water_u, skeleton_u) = column_base("drained"), column_base("undrained_a")
check(len(water_d) == 0, "drained ground reports no excess pore pressure")
check(len(water_u) == len(skeleton_u) and len(water_u) > 0,
      "undrained ground reports its excess pore pressure at every node (a Plastic phase reported "
      "none)")
if water_u:
    worst = max(abs(w + s - Q) for w, s in zip(water_u, skeleton_u))
    print(f"      undrained: the water takes {min(water_u):.3f} .. {max(water_u):.3f} kPa of the "
          f"{Q:.0f}, the skeleton {min(skeleton_u):.3f} .. {max(skeleton_u):.3f}")
    check(worst < 1e-6 * Q, "at every node the two add up to the surcharge")
    check(min(water_u) > 0.9 * Q, "and the water takes nearly all of it")
print(f"      base reaction added by the surcharge: drained {rise_d:.3f}, undrained {rise_u:.3f} kN/m "
      f"(q W = {Q * W:.1f})")
check(abs(rise_d / (Q * W) - 1) < 1e-6, "drained, the base takes the whole surcharge")
check(abs(rise_u / (Q * W) - 1) < 1e-6,
      "undrained, it takes the whole surcharge as well (it reported 3% of it while the reaction "
      "left the pore pressure out)")

# --- 2. the strip footing on undrained clay ----------------------------------------------------
SU, HALF = 20.0, 1.0
PRANDTL = (2.0 + math.pi) * SU * HALF


def footing(drainage):
    prj = katai.Project("footing", element="tri15", mesh_size=0.5, auto_refine=True)
    clay = prj.materials.mohr_coulomb("Clay", E=1.0e4, nu=0.3, c=SU, phi=0.0, psi=0.0, gamma=0.0,
                                      gamma_sat=0.0, drainage=drainage, tension_cutoff=False)
    prj.geometry.polygon([(0, 0), (6, 0), (6, 4), (0, 4)], material=clay,
                         fix=["full", "horizontal", "free", "horizontal"])
    d = prj.displacements.line((0, 4), (HALF, 4), uy=-0.05, name="Footing")
    prj.initial(procedure="k0", exclude=[d])
    prj.phases.plastic("push", activate=[d])
    job = katai.Job(prj.build())
    job.run()
    r = job.results()[-1]
    return -sum(r.reaction[2 * i + 1] for i in range(len(r.node_x))
                if abs(r.node_y[i] - 4) < 1e-9 and r.node_x[i] <= HALF + 1e-9)


f_d, f_u = footing("drained"), footing("undrained_a")
print(f"      footing force at 5 cm: drained {f_d:.2f}, undrained {f_u:.2f} kN/m "
      f"(Prandtl {PRANDTL:.2f})")
check(abs(f_u / PRANDTL - 1) < 0.04,
      f"the footing on undrained clay reads Prandtl's force within 4%: {f_u:.2f} against "
      f"{PRANDTL:.2f} (it read 22.7 from the effective stress alone)")
check(0.9 * PRANDTL < f_d < 1.02 * PRANDTL,
      f"the same clay declared drained is on its way to the same limit: {f_d:.2f}")

# --- 3. the joint under a block loaded drained and undrained -----------------------------------
CI, PHII = 3.5, 18.0
TANI = math.tan(math.radians(PHII))
PUSH = 0.01
T_NONE = (CI + GAMMA * H * TANI) * W
T_LOADED = (CI + (GAMMA * H + Q) * TANI) * W


def base_friction(drainage, surcharge):
    prj = katai.Project("shear box", element="tri6", mesh_size=0.25, auto_refine=False)
    block = prj.materials.linear_elastic("Block", E=2.0e4, nu=0.3, gamma=GAMMA, gamma_sat=GAMMA,
                                         drainage=drainage)
    rough = prj.materials.mohr_coulomb("Joint", E=2.0e4, nu=0.3, c=CI, phi=PHII, psi=0.0,
                                       gamma=GAMMA, gamma_sat=GAMMA, tension_cutoff=False)
    prj.geometry.polygon([(0, 0), (W, 0), (W, H), (0, H)], material=block,
                         fix=["full", "free", "free", "free"])
    prj.structures.interface((0, 0), (W, 0), name="Joint", iface_material=rough)
    hold = [prj.displacements.line((0, 0), (0, H), ux=0.0, name="hold left"),
            prj.displacements.line((W, 0), (W, H), ux=0.0, name="hold right")]
    push = [prj.displacements.line((0, 0), (0, H), ux=PUSH, name="push left"),
            prj.displacements.line((W, 0), (W, H), ux=PUSH, name="push right")]
    q = prj.loads.line_load((0, H), (W, H), qy=-Q, name="q")
    prj.initial(procedure="k0", exclude=[q] + hold + push)
    prj.phases.plastic("confine", activate=hold)
    if surcharge:
        prj.phases.plastic("surcharge", activate=[q])
    prj.phases.plastic("push", activate=push, deactivate=hold)
    job = katai.Job(prj.build())
    job.run()
    r = job.results()[-1]
    # the force on the pushed sides is the friction under the block (nothing else holds it)
    return sum(r.reaction[2 * i] for i in range(len(r.node_x))
               if (abs(r.node_x[i]) < 1e-9 or abs(r.node_x[i] - W) < 1e-9)
               and abs(r.displacement[2 * i] - PUSH) < 1e-6 * PUSH)


t_none = base_friction("drained", False)
t_drained = base_friction("drained", True)
t_undrained = base_friction("undrained_a", True)
print(f"      friction under the block: no load {t_none:.2f} | loaded drained {t_drained:.2f} | "
      f"loaded undrained {t_undrained:.2f} kN/m     closed form {T_NONE:.2f} | {T_LOADED:.2f} | "
      f"{T_NONE:.2f} and a little")
check(abs(t_none / T_NONE - 1) < 0.04,
      f"with no load the block slides at c_i W + W_block tan(phi_i): {t_none:.2f} against {T_NONE:.2f}")
check(abs(t_drained / T_LOADED - 1) < 0.04,
      f"a drained surcharge adds q W tan(phi_i): {t_drained:.2f} against {T_LOADED:.2f}")
check(1.0 <= t_undrained / t_none < 1.10,
      f"an undrained surcharge adds next to nothing: {t_undrained:.2f} against {t_none:.2f} "
      f"(friction on the transmitted stress would give about {t_drained:.0f})")

sys.exit(1 if failures else 0)
