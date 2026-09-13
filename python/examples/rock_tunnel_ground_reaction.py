"""A circular tunnel in a Hoek-Brown rock mass, against the closed-form solution.

THE PROBLEM. A tunnel of radius R is driven in a rock mass under a hydrostatic in-situ
stress sigma_0. The excavation is modelled the way the convergence-confinement method
means it: the internal pressure p_i is let down from sigma_0 to zero in stages, and the
wall displacement at each stage traces the GROUND REACTION CURVE. Below a critical
pressure a plastic annulus opens around the opening and grows as the support pressure
falls.

Rock mass and stress state are case D1 of Table A1.1 in Hoek, Carranza-Torres, Diederichs
& Corkum (2008), "Integration of geotechnical and structural design in tunnelling", 56th
Annual Geotechnical Engineering Conference, University of Minnesota -- a graphitic
phyllite at 1100 m depth, from the Yacambu-Quibor tunnel in Venezuela.

THE ORACLE is derived here rather than quoted, so that nothing rests on a remembered
formula. Compression positive, radial equilibrium of an axisymmetric ring is

    d(sigma_r)/dr = (sigma_theta - sigma_r) / r

and inside the plastic annulus the stress pair sits ON the Hoek-Brown envelope,

    sigma_theta = sigma_r + sigma_ci (m_b sigma_r/sigma_ci + s)^a

Substituting the second into the first separates the variables and integrates in closed
form for any a != 1:

    ln(r/R) = 1/(m_b (1-a)) [ (m_b sigma_r/sigma_ci + s)^(1-a)
                              - (m_b p_i/sigma_ci + s)^(1-a) ]

At the elastic-plastic boundary the elastic field around a hole in an infinite medium
under hydrostatic sigma_0 gives sigma_r + sigma_theta = 2 sigma_0, so the radial stress
there solves

    2 (sigma_0 - sigma_r^ep) = sigma_ci (m_b sigma_r^ep/sigma_ci + s)^a

and putting that into the integral gives the plastic radius. Setting a = 1/2 reduces both
to the scaled forms of Carranza-Torres & Fairhurst (1999) -- R_pl = R exp(2(sqrt(S^ep) -
sqrt(P_i))) and P_i^cr = (1/16)[1 - sqrt(1 + 16 S_0)]^2 -- which is how the derivation is
checked against the literature rather than against itself.

THE MODEL IS A RADIAL STRIP IN AXISYMMETRY, and that is worth a sentence because it looks
like a different problem. In axisymmetry the three normal stresses are (sigma_r, sigma_z,
sigma_theta). Map the tunnel onto it: r is the radius, theta the hoop, and z the TUNNEL
AXIS. Holding u_z = 0 on both horizontal edges is then exactly the tunnel's plane-strain
condition, and the whole cross-section reduces to a one-dimensional strip with straight
edges -- no polygonal approximation of a circle, and two pressure boundaries instead of
forty facet loads.

Units: kN, m. Run with the built package on PYTHONPATH::

    python rock_tunnel_ground_reaction.py
"""
import math
import sys

import katai

# ---- the case, from Table A1.1 (case D1) --------------------------------------
SIGMA_CI = 50.0e3      # kPa, uni-axial compressive strength of the INTACT rock
MI = 7.0               # intact rock parameter
GSI = 48.0             # Geological Strength Index
D = 0.0                # disturbance factor
E_RM = 7.5e6           # kPa, rock mass Young's modulus
NU = 0.25
SIGMA_0 = 28.0e3       # kPa, hydrostatic in-situ stress
R = 2.5                # m, tunnel radius
B = 100.0              # m, outer radius of the model (40 R: the far field is far)
H = 2.0                # m, strip height -- only the plane-strain constraint needs it

# Support pressures the excavation is let down through [kPa]. Coarse while the rock is
# still elastic, close together once the plastic annulus starts to open, because that is
# where the answer changes fastest -- and where the equilibrium iteration works hardest.
LEVELS = [24500.0, 21000.0, 17500.0, 15000.0, 14000.0, 13000.0, 12000.0, 10000.0,
          8000.0, 6000.0, 4000.0, 3000.0, 2000.0, 1000.0, 500.0, 250.0,
          150.0, 100.0, 60.0, 30.0, 15.0, 5.0, 0.0]


# ---- the closed form ----------------------------------------------------------
def rock_mass_constants():
    """m_b, s and a from GSI and D (Hoek, Carranza-Torres & Corkum 2002)."""
    mb = MI * math.exp((GSI - 100.0) / (28.0 - 14.0 * D))
    s = math.exp((GSI - 100.0) / (9.0 - 3.0 * D))
    a = 0.5 + (math.exp(-GSI / 15.0) - math.exp(-20.0 / 3.0)) / 6.0
    return mb, s, a


def critical_pressure():
    """Radial stress at the elastic-plastic boundary; also the internal pressure below
    which a plastic annulus exists at all (there R_pl = R)."""
    mb, s, a = rock_mass_constants()

    def g(x):
        return 2.0 * (SIGMA_0 - x) - SIGMA_CI * (mb * x / SIGMA_CI + s) ** a

    lo, hi = 0.0, SIGMA_0
    if g(lo) < 0.0:
        return None                      # the rock never yields, even unsupported
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if g(mid) > 0.0 else (lo, mid)
    return 0.5 * (lo + hi)


def plastic_radius(p_i):
    mb, s, a = rock_mass_constants()
    x_ep = critical_pressure()
    if x_ep is None or p_i >= x_ep:
        return R
    k = 1.0 / (mb * (1.0 - a))
    return R * math.exp(k * ((mb * x_ep / SIGMA_CI + s) ** (1.0 - a)
                             - (mb * p_i / SIGMA_CI + s) ** (1.0 - a)))


def radial_stress(r, p_i):
    """Closed-form sigma_r at radius r, compression positive."""
    mb, s, a = rock_mass_constants()
    r_pl, x_ep = plastic_radius(p_i), critical_pressure()
    if r >= r_pl:                                    # elastic, infinite medium
        return SIGMA_0 - (SIGMA_0 - x_ep) * (r_pl / r) ** 2
    k = 1.0 / (mb * (1.0 - a))                       # plastic branch, inverted
    v = (mb * p_i / SIGMA_CI + s) ** (1.0 - a) + math.log(r / R) / k
    return (v ** (1.0 / (1.0 - a)) - s) * SIGMA_CI / mb


# ---- the finite element model -------------------------------------------------
# The settings the published comparison was measured at. `--quick` is the same problem on a
# coarser mesh and a shorter list of support pressures: it is what the test suite runs, because
# an example nobody runs stops working quietly, and a five-minute one does not get run.
QUICK = "--quick" in sys.argv
MESH = 1.0 if QUICK else 0.4
STEPS = 40 if QUICK else 90
if QUICK:
    LEVELS = [20000.0, 14000.0, 12000.0, 10000.0, 8000.0, 6000.0, 4500.0, 3000.0]


def build(mesh_size=MESH, steps=STEPS):
    prj = katai.Project("Rock tunnel, ground reaction curve", kind="axisymmetric",
                        mesh_size=mesh_size, auto_refine=False)
    # The dilatancy is left at zero: the closed form above is a stress solution and does
    # not read the flow rule, but the WALL DISPLACEMENT does, and zero is the usual
    # assumption for a rock mass at this strength.
    rock = prj.materials.hoek_brown("Graphitic phyllite", E=E_RM, nu=NU,
                                    sigma_ci=SIGMA_CI, mi=MI, gsi=GSI, D=D,
                                    psi=0.0, sigma_psi=0.0, gamma=0.001)
    # r = R is the tunnel wall, r = B the far field; u_z = 0 top and bottom is plane
    # strain along the tunnel axis.
    prj.geometry.rectangle(R, 0.0, B, H, material=rock,
                           bottom="vertical", top="vertical",
                           left="free", right="free")
    prj.loads.line_load((B, 0.0), (B, H), qx=-SIGMA_0, name="In-situ stress")
    # The support pressure that never comes off. The relief pieces below sum to
    # sigma_0 - LEVELS[-1], so without this the tunnel would start at that reduced pressure and
    # end at ZERO rather than at the last level -- the stages would be right and every pressure
    # they are compared against would be wrong. It is zero for the published list, which ends at
    # zero, and it is what makes the shortened `--quick` list mean what it says.
    if LEVELS[-1] > 0.0:
        prj.loads.line_load((R, 0.0), (R, H), qx=+LEVELS[-1], name="Support that stays")
    # The rest as one load per stage, so that letting it down is a plain deactivation and each
    # phase's displacement is that stage's convergence.
    pieces, previous = [], SIGMA_0
    for target in LEVELS:
        pieces.append(prj.loads.line_load((R, 0.0), (R, H), qx=+(previous - target),
                                          name=f"Relief to {target / 1000:.2f} MPa"))
        previous = target
    for piece, target in zip(pieces, LEVELS):
        prj.phases.plastic(f"p_i = {target / 1000:.2f} MPa", deactivate=[piece],
                           steps=steps, max_iterations=1500)
    return prj


def wall_displacement(result):
    xs, ys, u = result.node_x, result.node_y, result.displacement
    best, index = 1e30, -1
    for i in range(len(xs)):
        d = abs(xs[i] - R) + abs(ys[i] - 0.5 * H)
        if d < best:
            best, index = d, i
    return u[2 * index]


def main():
    mb, s, a = rock_mass_constants()
    print(f"rock mass:  m_b = {mb:.4f}   s = {s:.6f}   a = {a:.4f}")
    print(f"            (Table A1.1 of the source prints 1.093, 0.0031, 0.507)")
    print(f"critical support pressure = {critical_pressure() / 1000:.3f} MPa")
    print(f"plastic radius when unsupported = {plastic_radius(0.0) / R:.3f} R\n")

    job = katai.Job(build().build())
    ok = job.run()
    results = job.results()

    print("   p_i [MPa]   wall closure [mm]   R_pl/R (closed form)")
    total = 0.0
    for k, p_i in enumerate(LEVELS):
        if k + 1 >= len(results) or not results[k + 1].ok:
            print(f"   -- the run stops here; see the phase message")
            break
        total += wall_displacement(results[k + 1])
        print(f"   {p_i / 1000:9.2f}   {abs(total) * 1000:17.2f}   {plastic_radius(p_i) / R:19.3f}")

    last = next(r for r in reversed(results) if r.ok)
    reached = min(LEVELS[k] for k in range(len(LEVELS))
                  if k + 1 < len(results) and results[k + 1].ok)
    xs, ys, st = last.node_x, last.node_y, last.stress
    profile = sorted((xs[i], -st[i][0]) for i in range(len(xs))
                     if abs(ys[i] - 0.5 * H) < 1e-6)
    print(f"\n   radial stress at p_i = {reached / 1000:.3f} MPa"
          f"   (closed-form R_pl = {plastic_radius(reached) / R:.3f} R)")
    print("      r/R    FE [MPa]   closed form [MPa]      difference")
    for radius, sigma_r in profile:
        if radius > 6.0 * R:
            continue
        exact = radial_stress(radius, reached)
        print(f"   {radius / R:6.2f}   {sigma_r / 1000:9.3f}   {exact / 1000:16.3f}"
              f"   {100.0 * (sigma_r - exact) / max(abs(exact), 1.0):+13.2f}%")
    if not ok:
        print("\n   (the last stage did not converge -- the phase message says why)")

    # ---- what makes this an example that cannot rot quietly --------------------
    # Three claims, each of which fails the RUN rather than the reader.
    failures = []
    # 1. While the rock is elastic every stage closes the wall by (1+nu) dp R / E exactly.
    elastic = [k for k, q in enumerate(LEVELS) if plastic_radius(q) <= R * (1.0 + 1e-9)]
    for k in elastic[:4]:
        if k + 1 >= len(results) or not results[k + 1].ok:
            break
        want = (1.0 + NU) * ((LEVELS[k - 1] if k else SIGMA_0) - LEVELS[k]) * R / E_RM
        got = abs(wall_displacement(results[k + 1]))
        if abs(got / want - 1.0) > 0.02:
            failures.append(f"elastic stage {k + 1}: {got*1000:.3f} mm, expected "
                            f"{want*1000:.3f} mm")
    # 2. The radial stress lands on the closed form on both sides of the boundary.
    inside = outside = 0
    for radius, sigma_r in profile:
        if radius <= R + 1e-9 or radius > 6.0 * R:
            continue                    # the wall node is reported, never asserted
        exact = radial_stress(radius, reached)
        rel = abs(sigma_r - exact) / max(abs(exact), 1.0)
        if radius < plastic_radius(reached):
            inside += 1
        else:
            outside += 1
        if rel > 0.02:
            failures.append(f"sigma_r at r/R = {radius/R:.2f}: {rel*100:.2f}% from the "
                            f"closed form")
    # 3. ...and it actually SAMPLED both sides. A band on a region with no samples in it is
    #    an assertion that cannot fail, which is worse than no assertion at all.
    if inside < 1 or outside < 3:
        failures.append(f"the profile sampled {inside} node(s) inside the plastic annulus "
                        f"and {outside} outside: too few to be asserting anything")
    if failures:
        print("\nFAILED:")
        for f in failures:
            print("   " + f)
        raise SystemExit(1)
    print(f"\nOK: {len(elastic)} elastic stage(s) match (1+nu) dp R / E, and the radial "
          f"stress is the closed form's at {inside + outside} nodes on both sides")


if __name__ == "__main__":
    main()
