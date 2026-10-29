"""The water at a wall's joint: the joint's law is a law of the soil skeleton, the pore water goes past it.

A wall with interfaces on both sides, wished in place in ground whose water table is at the surface:
  (1) installing it moves nothing -- the ground stood in equilibrium before, and the wall takes over
      the stress that was there, effective stress and pore water alike;
  (2) the joint holds the EFFECTIVE normal stress: the same as in the same ground dry, with the
      buoyant unit weight gamma' for its weight -- the ground of the same effective stresses;
  (3) it rubs with that stress: pulled out of the ground, it holds 0.9 of the friction its effective
      stress gives along both faces, c_i - sigma'_n tan(phi_i) summed over the joint, and not 1.1.
      With the pore pressure in the friction it held three times as much;
  (4) a pit dug in front of it below the water table, the water left standing in the pit: the wall
      bends exactly as in the dry ground of buoyant weight -- the water against its two faces
      cancels, and the pit floor's water bears on the floor alone.
Until 2026-10-29 the joint took the pore pressure into its own normal stress: the installation moved
the ground by p / k_n at every wall below the water table and the friction read the total stress.
"""
import sys

import katai

G = 9.81
failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def model(wet, *, pull=0.0, pit=False):
    prj = katai.Project("joint water", element="tri6", mesh_size=1.0, auto_refine=False)
    gam = 18.0 if wet else 18.0 - G
    soil = prj.materials.mohr_coulomb("Soil", E=2.0e4, nu=0.3, c=1.0, phi=25.0, psi=0.0,
                                      gamma=gam, gamma_sat=18.0 if wet else gam)
    M = prj._materials[soil.index]
    M.rinter_rigid = False
    M.Rinter = 0.5
    F, H, X = "free", "horizontal", "full"
    prj.geometry.polygon([(0, -20), (20, -20), (20, -4), (20, 0), (0, 0)], material=soil, name="Behind",
                         fix=[X, F, F, F, H])
    dig = prj.geometry.polygon([(20, -4), (40, -4), (40, 0), (20, 0)], material=soil, name="Pit",
                               fix=[F, H, F, F])
    prj.geometry.polygon([(20, -20), (40, -20), (40, -4), (20, -4)], material=soil, name="Below",
                         fix=[X, H, F, F])
    wall = prj.structures.plate((20, 0), (20, -10), EA=1.0e7, EI=1.0e5, w=0.0, nu=0.0,
                                interfaces="both", name="Wall")
    up = prj.loads.point_load((20, 0), qy=pull, name="Pull") if pull else None
    prj.water.table(0.0 if wet else -100.0)
    prj.initial(procedure="k0", exclude=[wall] + ([up] if up else []))
    prj.phases.plastic("Install", activate=[wall])
    if up:
        prj.phases.plastic("Pull", activate=[up])
    if pit:
        prj.phases.plastic("Dig", deactivate=[dig])
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job.results()


def joint_stations(r):
    out = []
    for itf in r.interface_forces:
        out += [(q.x, q.y, q.sigma_n) for q in itf.stations]
    return sorted(out)


# (1) + (2): install
ok_w, res_w = model(True)
ok_d, res_d = model(False)
r_w, r_d = res_w[1], res_d[1]
check(ok_w and ok_d, "both installations run")
check(r_w.max_disp < 1e-9, f"installing the wall below the water table moves nothing "
                           f"(max |u| = {r_w.max_disp:.2e} m)")
sw, sd = joint_stations(r_w), joint_stations(r_d)
dev = max(abs(a[2] - b[2]) for a, b in zip(sw, sd)) if len(sw) == len(sd) and sw else 1e9
big = max(abs(b[2]) for b in sd) if sd else 0.0
print(f"  (2) {len(sw)} joint stations, max |sigma_n(wet) - sigma_n(dry, gamma')| = {dev:.2e} of {big:.1f} kPa")
check(len(sw) > 0 and dev < 1e-6 * big, "the joint holds the effective normal stress, as in the "
                                        "dry ground of buoyant weight")

# (3): pull-out against the joint's own effective normal stress
tan_i = 0.5 * 0.46630765815499864   # R_inter tan(25 deg)
c_i = 0.5 * 1.0
cap = 0.0
for itf in r_w.interface_forces:
    st = sorted(itf.stations, key=lambda q: q.s)
    for a, b in zip(st[:-1], st[1:]):
        ta = max(0.0, c_i - a.sigma_n * tan_i)
        tb = max(0.0, c_i - b.sigma_n * tan_i)
        cap += 0.5 * (ta + tb) * (b.s - a.s)
print(f"  (3) friction of the joint on its effective normal stress, both faces: {cap:.2f} kN/m")
# Pulled with more than that the joint slips along its length (the ground under the toe, which the
# joint does not reach, takes the rest), and where it slips it holds c_i - sigma'_n tan(phi_i).
ok_p, res_p = model(True, pull=1.5 * cap)
slipping, worst = 0, 0.0
stations = [q for itf in res_p[-1].interface_forces for q in itf.stations if -9.5 < q.y < -0.5]
for q in stations:
    if q.slipping:
        slipping += 1
        lim = c_i - q.sigma_n * tan_i
        worst = max(worst, abs(abs(q.tau) - lim) / lim)
print(f"  (3) pulled with 1.5 x that: {slipping} of {len(stations)} joint points slipping, "
      f"worst |tau| against c_i - sigma'_n tan(phi_i): {worst:.2e}")
check(ok_p and slipping >= 0.9 * len(stations) and worst < 1e-3,
      "pulled out, the joint slips with the friction of its effective normal stress")

# (4): a submerged pit against the dry ground of buoyant weight
ok_pw, res_pw = model(True, pit=True)
ok_pd, res_pd = model(False, pit=True)
mw = max(abs(s.M) for s in res_pw[-1].struct_forces[0].stations)
md = max(abs(s.M) for s in res_pd[-1].struct_forces[0].stations)
uw = max(abs(s.ux) for s in res_pw[-1].struct_forces[0].stations)
ud = max(abs(s.ux) for s in res_pd[-1].struct_forces[0].stations)
print(f"  (4) wall max |M|: submerged pit {mw:.4f}, dry ground of gamma' {md:.4f} kNm/m; "
      f"max |ux| {1000 * uw:.4f} / {1000 * ud:.4f} mm")
check(ok_pw and ok_pd and md > 1.0 and abs(mw - md) < 1e-3 * md and abs(uw - ud) < 1e-3 * ud,
      "a pit dug below the water table bends the wall as the dry ground of buoyant weight does (0.1%)")

sys.exit(1 if failures else 0)
