"""Pore pressures of a staged phase from a steady groundwater flow (Phase.pore_from_flow).

Three checks, each against a closed form or against an independent path through the program:
  (1) with one water table and no wall the flow is hydrostatic, so the phase takes the same pore
      pressures the phreatic line gives -- to round-off -- and moves next to nothing;
  (2) a full cut-off wall between two water tables, impermeable down to the closed base: each side
      is hydrostatic at its own table, everywhere, and the wall carries the same bending moment as the
      same wall under the phreatic line with the same two tables -- the difference of water pressure
      on its two faces is in its forces, which a head handed over from an unsplit mesh could not do;
  (3) the same wall stopping half way: water flows under its toe, so in front of it the pore
      pressure is above the low side's hydrostatic value and behind it below the high side's;
  (4) a pit dug in front of that wall and pumped dry to its floor: water rises through the floor,
      and the phase reports the hydraulic-heave ratio (K2D-A024): the least over every column and
      depth, found in the column beside the wall at the depth of its toe -- the classical one,
      gamma' / (gamma_w i_avg) read by hand from the heads on the wall's face over its embedment,
      to 5% (the first column stands half an element off the wall: 4.549 against 4.558 at 1 m).
      Before the wall's floor node took its own side's table this was 2.13 and fell under 1 as
      the mesh was refined.
"""
import re
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


def model(wall_bottom, flow, left=-2.0, right=-8.0):
    prj = katai.Project("flow pore", element="tri6", mesh_size=1.0, auto_refine=False)
    m = prj.materials.linear_elastic("Soil", E=5e4, nu=0.3, gamma=18.0, gamma_sat=20.0, kx=1.0, ky=1.0)
    F, H, X = "free", "horizontal", "full"
    prj.geometry.polygon([(0, -20), (15, -20), (15, 0), (0, 0)], material=m, name="L", fix=[X, F, F, H])
    prj.geometry.polygon([(15, -20), (30, -20), (30, 0), (15, 0)], material=m, name="R", fix=[X, H, F, F])
    if wall_bottom is not None:
        prj.structures.plate((15, 0), (15, wall_bottom), EA=1e7, EI=1e5, w=0.0, interfaces="both",
                             barrier="impermeable", name="Wall")
    prj.water.table(left)
    prj.initial(procedure="k0")
    table = left if right == left else [(0, left), (14.999, left), (15.001, right), (30, right)]
    prj.phases.plastic("Water", water=table, pore_pressure="flow" if flow else "phreatic")
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job.results()[-1]


def heads(r):
    x, y, p = np.array(r.node_x), np.array(r.node_y), np.array(r.pore)
    return x, y, p, y + p / G


# (1) one table, no wall
ok, r = model(None, True, -2.0, -2.0)
x, y, p, h = heads(r)
dev = float(np.abs(p - G * np.maximum(0.0, -2.0 - y)).max())
print(f"  (1) max |u - u_hydrostatic| = {dev:.2e} kPa, phase max|u| = {1000 * r.max_disp:.4f} mm")
check(ok and dev < 1e-9, "one water table and no wall: the flow is hydrostatic")
check(r.max_disp < 1e-4, "and the phase moves next to nothing (the two assemblies differ only where "
                         "elements straddle the table)")

# (2) full cut-off
ok, rf = model(-20.0, True)
x, y, p, h = heads(rf)
sat = p > 1e-6
left_dev = float(np.abs(h[sat & (x < 14.9)] + 2.0).max())
right_dev = float(np.abs(h[sat & (x > 15.1)] + 8.0).max())
print(f"  (2) head off its table: behind {left_dev:.2e} m, in front {right_dev:.2e} m")
check(ok and left_dev < 1e-6 and right_dev < 1e-6,
      "a full cut-off to the closed base: each side hydrostatic at its own table")
okp, rp = model(-20.0, False)
mf = max(abs(s.M) for s in rf.struct_forces[0].stations)
mp = max(abs(s.M) for s in rp.struct_forces[0].stations)
print(f"  (2) wall max |M|: flow {mf:.3f}, phreatic {mp:.3f} kNm/m")
check(okp and mp > 1.0 and abs(mf - mp) < 0.01 * mp,
      "and the wall carries the moment of the same two tables under the phreatic line (1%)")

# (3) partial cut-off
ok, r3 = model(-10.0, True)
x, y, p, h = heads(r3)
front = (x > 15.1) & (x < 17.0) & (y < -11.0) & (y > -14.0)
behind = (x < 14.9) & (x > 13.0) & (y < -11.0) & (y > -14.0)
hf, hb = float(h[front].mean()), float(h[behind].mean())
print(f"  (3) mean head under the toe level: in front {hf:.3f} m, behind {hb:.3f} m")
check(ok and -8.0 + 0.1 < hf < -2.0 and -8.0 < hb < -2.0 - 0.1 and hb > hf,
      "a wall stopping half way: water flows under its toe, raising the head in front of it and "
      "drawing it down behind")

# (4) a pit pumped dry in front of a wall stopping half way
def pit():
    prj = katai.Project("flow heave", element="tri6", mesh_size=1.0, auto_refine=False)
    m = prj.materials.linear_elastic("Soil", E=5e4, nu=0.3, gamma=18.0, gamma_sat=20.0, kx=1.0, ky=1.0)
    F, H, X = "free", "horizontal", "full"
    prj.geometry.polygon([(0, -20), (15, -20), (15, 0), (0, 0)], material=m, name="L", fix=[X, F, F, H])
    prj.geometry.polygon([(15, -20), (30, -20), (30, -4), (15, -4)], material=m, name="R",
                         fix=[X, H, F, F])
    dig = prj.geometry.polygon([(15, -4), (30, -4), (30, 0), (15, 0)], material=m, name="Pit",
                               fix=[F, H, F, F])
    prj.structures.plate((15, 0), (15, -8), EA=1e7, EI=1e5, w=0.0, interfaces="both",
                         barrier="impermeable", name="Wall")
    prj.water.table(-1.0)
    prj.initial(procedure="k0")
    prj.phases.plastic("Dig", deactivate=[dig], dry=[dig],
                       water=[(0, -1.0), (14.999, -1.0), (15.001, -4.0), (30, -4.0)],
                       pore_pressure="flow")
    job = katai.Job(prj.build())
    ok = job.run()
    return ok, job.results()[-1]


ok, r4 = pit()
x, y, p, h = heads(r4)
notes = [d for d in r4.diagnostics if d.code in ("K2D-A024", "K2D-A025")]
fac = None
if notes:
    mm = re.search(r"(?:is|weight over excess water pressure) (-?[0-9][0-9.eE+-]*)", notes[0].message)
    fac = float(mm.group(1)) if mm else None


def head_near(px, py):
    # on the pit's side of the wall only: the wall's twin nodes share its coordinates
    d = np.where(x > 15.0 + 1e-6, (x - px) ** 2 + (y - py) ** 2, np.inf)
    return h[int(np.argmin(d))]


i_avg = (head_near(15.0, -8.0) - head_near(15.0, -4.0)) / 4.0
f_hand = (20.0 - G) / (G * i_avg)
print(f"  (4) {notes[0].code if notes else 'no note'}: ratio {fac}; by hand next to the wall over its "
      f"embedment: i_avg {i_avg:.4f}, gamma'/(gamma_w i_avg) = {f_hand:.4f}")
check(ok and fac is not None and abs(fac - f_hand) <= 0.05 * f_hand,
      "a pit pumped dry in front of the wall: the heave ratio is reported, and is the classical one "
      "of the column next to the wall over its embedment (5%)")

sys.exit(1 if failures else 0)
