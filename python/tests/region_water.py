"""Region water (Phase.poly_water): a soil region's pore pressure from a level of its own, or
interpolated with depth between the regions above and below it.

Checked against closed forms on a 10 m wide column (y = -12..0) whose top 2 m are dug out:
  (1) the region between the dug floor (-2 m, the pit pumped to its floor) and -6 m interpolated:
      u = u(-6) (y_top - y)/(y_top - y_bot) with u(-6) = gamma_w (-1 - (-6)) from the water table
      at -1 m below it, so the pore pressure rises linearly from 0 at the floor to the table's
      value at -6, and below -6 it is the table's hydrostatic value; every node to 1e-4 kPa (each
      node is read a millionth of the way into the element it is read for, so a node on a region's
      boundary falls inside the region: about 1e-6 m, 1e-5 kPa);
  (2) the same region with a level of its own at -3 m: hydrostatic from -3 m inside it (1e-4 kPa);
  (3) the weight of the dug soil leaves the floor whichever way the pore pressure is given: with
      the water table deep below the dig, the floor heaves the one-dimensional unloading
      gamma h H / E_oed whether the phase reads the water table or a region rule (a head field) --
      until 2026-10 the head-field path kept part of the dug soil's weight on the floor (10.01
      against 10.70 mm);
  (4) a region still on the water table is the water table: the phase is the phreatic one.
"""
import sys

import numpy as np

import katai

G = 9.81
failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def column(table, region_water=None, E=2.0e4, nu=0.3):
    prj = katai.Project("region water", element="tri6", mesh_size=1.0, auto_refine=False)
    m = prj.materials.linear_elastic("S", E=E, nu=nu, gamma=18.0, gamma_sat=20.0, kx=1.0, ky=1.0)
    pit = prj.geometry.polygon([(0, -2), (10, -2), (10, 0), (0, 0)], material=m, name="Pit")
    mid = prj.geometry.polygon([(0, -6), (10, -6), (10, -2), (0, -2)], material=m, name="Mid")
    prj.geometry.polygon([(0, -12), (10, -12), (10, -6), (0, -6)], material=m, name="Base")
    prj.water.table(table)
    prj.initial(procedure="k0")
    rw = None if region_water is None else {mid: region_water}
    prj.phases.plastic("Dig", deactivate=[pit], dry=[pit], region_water=rw)
    job = katai.Job(prj.build())
    ok = job.run()
    r = job.results()[-1]
    return ok, r, np.array(r.node_x), np.array(r.node_y), np.array(r.pore)


# (1) interpolated between the pumped floor and the table below
ok, r, x, y, p = column(-1.0, "interpolate")
mid_n = (y < -2.0 - 1e-9) & (y > -6.0 + 1e-9)
base_n = y < -6.0 - 1e-9
u6 = G * 5.0
exact_mid = u6 * (-2.0 - y) / 4.0
dev_mid = float(np.abs(p[mid_n] - exact_mid[mid_n]).max())
dev_base = float(np.abs(p[base_n] - G * (-1.0 - y[base_n])).max())
print(f"  (1) interpolated: max |u - closed form| {dev_mid:.2e} kPa in it, {dev_base:.2e} kPa below")
check(ok and dev_mid < 1e-4 and dev_base < 1e-4 and mid_n.sum() > 10,
      "interpolated between the pumped floor and the water table below: linear (1e-4 kPa)")
check(any(d.code == "K2D-A024" for d in r.diagnostics),
      "and the heave check reads the field it states (water rising to the floor)")

# (2) a level of its own
ok, r, x, y, p = column(-1.0, -3.0)
dev = float(np.abs(p[mid_n] - G * np.maximum(0.0, -3.0 - y[mid_n])).max())
print(f"  (2) own level -3 m: max |u - closed form| {dev:.2e} kPa")
check(ok and dev < 1e-4, "a level of its own: hydrostatic from it (1e-4 kPa)")

# (3) the dug soil's weight leaves the floor on the head-field path too
E, nu = 2.0e4, 0.3
eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
exact = 18.0 * 2.0 * 10.0 / eoed * 1000.0


def floor_heave(rw):
    ok, r, x, y, p = column(-20.0, rw)
    d = np.array(r.displacement)
    k = int(np.argmin((x - 5.0) ** 2 + (y + 2.0) ** 2))
    return ok, d[2 * k + 1] * 1000.0


ok_t, h_t = floor_heave(None)
ok_r, h_r = floor_heave(-25.0)       # a region rule, the water still deep below the dig
print(f"  (3) floor heave: water table {h_t:.4f} mm, region rule {h_r:.4f} mm, 1D unloading {exact:.4f} mm")
check(ok_t and ok_r and abs(h_r - h_t) < 1e-6 * exact and abs(h_t - exact) < 0.01 * exact,
      "the dug soil's weight leaves the floor on the head-field path as on the phreatic one")

# (4) a region on the water table is the water table
ok_a, ra, xa, ya, pa = column(-1.0, None)
ok_b, rb, xb, yb, pb = column(-1.0, "table")
check(ok_a and ok_b and float(np.abs(pa - pb).max()) == 0.0,
      'region_water "table" is the water table itself, bit for bit')

sys.exit(1 if failures else 0)
