"""A coupled phase reads its force balance and its volume balance each against its own scale.

One norm over the two, with a floor tied to the force the soil already carries, let a time step close at its first
iterate once the flow of the step fell under that floor, with nothing moved. The pressure and the settlement froze and
the phase reported itself solved: 18% short of the end settlement on a column of Mohr-Coulomb sand seated at 1000 kPa,
the excess pressure stopped at 2.9 kPa where it should be gone, and more of it on more steps. Against Terzaghi's column, on
the same model through the consolidation solver as well. The change of permeability is read by both.
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


H, K, GW, E, NU, DQ = 2.0, 1e-3, 9.81, 2.0e4, 0.3, 10.0
E_INIT = 0.5
KW = 2.0e6
EOED = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
MV = 1.0 / EOED
N_OVER_KW = (E_INIT / (1 + E_INIT)) / KW


def terzaghi_u(tv, nterm=200):
    return 1.0 - sum(2.0 / m ** 2 * math.exp(-m ** 2 * tv) for m in ((2 * j + 1) * math.pi / 2.0 for j in range(nterm)))


def column(route, seat, steps, tv_end=2.0, dq=DQ, c=1.0e4, ck=0.0, e=E, nu=NU, duration=None):
    """A weightless column, drained at its top, seated on `seat` kPa, with `dq` more at the first instant."""
    prj = katai.Project(f"column {route}", element="tri6", mesh_size=0.25, auto_refine=False, corner_refine=False)
    sand = prj.materials.mohr_coulomb("Sand", E=e, nu=nu, c=c, phi=30.0, psi=0.0, gamma=0.0, gamma_sat=0.0, k=K, e_init=E_INIT,
                                      ck=ck, k0=1.0, tension_cutoff=False)
    prj.geometry.rectangle(0, 0, 1, H, material=sand, bottom="full", right="horizontal", top="free", left="horizontal",
                           flow={"top": ("head", H)})
    prj.loads.line_load((0, H), (1, H), qy=-seat, name="seat")
    added = prj.loads.line_load((0, H), (1, H), qy=-dq, name="dq")
    prj.initial(procedure="gravity", exclude=[added])
    eoed = e * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    tend = duration if duration else tv_end * H ** 2 / (K * eoed / GW)
    if route == "consolidation":
        prj.phases.consolidation("c", duration=tend, steps=steps, activate=[added], load_at_start=True)
    else:
        prj.phases.fully_coupled("c", duration=tend, steps=steps, activate=[added])
    job = katai.Job(prj.build())
    ok = job.run()
    res = job.results()
    if not ok or not res:
        print("   (" + job.message()[:200] + ")")
        return None
    r = res[-1]
    return list(r.consol_time), list(r.consol_settlement), list(r.consol_excess_pore)


for route, seat, steps in (("fully_coupled", 100.0, 200), ("fully_coupled", 1000.0, 200), ("fully_coupled", 100.0, 800),
                           ("consolidation", 1000.0, 200)):
    got = column(route, seat, steps)
    check(got is not None, f"{route}, seated at {seat:.0f} kPa, {steps} steps: the phase closes")
    if got is None:
        continue
    t, s, p = got
    u0 = DQ * MV / (MV + N_OVER_KW)
    s0 = (DQ - u0) * H * MV
    cv = K / (GW * (MV + N_OVER_KW))
    exact = s0 + (DQ * H * MV - s0) * terzaghi_u(cv * t[-1] / H ** 2)
    err = 100.0 * (s[-1] / exact - 1.0)
    print(f"   {route:14s} seat {seat:6.0f} steps {steps:4d}: settlement {s[-1]:.6e} m ({err:+.3f}%), excess pressure {p[-1]:.3f} kPa")
    check(abs(err) < 0.2, f"{route}, seated at {seat:.0f} kPa, {steps} steps: the end settlement (0.2%)")
    check(p[-1] < 0.2, f"{route}, seated at {seat:.0f} kPa, {steps} steps: the pressure dissipates")

# The change of permeability: the same column, loaded enough to change its void ratio.
runs = {}
for route in ("consolidation", "fully_coupled"):
    for ck in (0.0, 0.05):
        eoed = 2000.0
        got = column(route, 100.0, 100, dq=50.0, ck=ck, e=2000.0, nu=0.0, duration=0.6 * H ** 2 / (K * eoed / GW))
        check(got is not None, f"{route}, ck = {ck}: the phase closes")
        if got:
            runs[(route, ck)] = got
            print(f"   {route:14s} ck {ck:4.2f}: settlement {got[1][-1]:.5f} m, excess pressure {got[2][-1]:.2f} kPa")
if len(runs) == 4:
    s_c0, s_c = runs[("consolidation", 0.0)][1][-1], runs[("consolidation", 0.05)][1][-1]
    s_f0, s_f = runs[("fully_coupled", 0.0)][1][-1], runs[("fully_coupled", 0.05)][1][-1]
    check(s_f < 0.5 * s_f0, "the coupled phase reads the change of permeability: the column is held back by it")
    check(abs(s_f / s_c - 1.0) < 0.05, "and lands with the consolidation solver (5%)")

print("\n" + ("OK" if not failures else f"{failures} check(s) failed"))
sys.exit(1 if failures else 0)
