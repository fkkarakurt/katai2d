"""Liquefaction triggering by the simplified method of TBDY 2018, Annex 16B.

A design-code evaluation beside the finite elements, not through them: SPT blow counts
corrected to N1,60 and for fines, the cyclic resistance of a magnitude 7.5 earthquake scaled
to the design magnitude, against the cyclic shear stress of the simplified method; a layer is
safe where tau_R / tau_eq >= 1.10 (Eq. 16.3).

    import katai.liquefaction as liq
    rows = liq.boring(tests=[(3.0, 16, 5), (6.0, 20, 16)], water_depth=1.8,
                      gamma=17.0, gamma_sat=20.0, SDS=0.75, Mw=7.5,
                      rod_factor=lambda z: 0.95)
    for r in rows: print(r["depth"], r["verdict"], r["FS"])

What the code does not say is not filled in: a test whose N1,60f reaches 30 is reported as
beyond the curve of Eq. 16B.4b, with no factor of safety.
"""
from . import _core

__all__ = ["spt", "boring", "rod_factor_tbdy"]


def rod_factor_tbdy(rod_length):
    """C_R from the rod length (Table 16B.1): 0.75 below 4 m, 0.85 to 6 m, 0.95 to 10 m, 1.00
    from 10 m (a boundary goes to the deeper row). The table starts at 3 m; a shorter rod takes
    0.75, its first row."""
    return _core.liquefaction_rod_factor(float(rod_length))


def spt(depth, N, fines, sigma_v, sigma_v_eff, *, SDS, Mw=7.5, CR=1.0, CS=1.0, CB=1.0, CE=1.0,
        below_water=True, susceptible=True, sigma_v_eff_at_test=0.0, required=1.10):
    """One SPT test. Stresses in kPa, depth in m, fines content in %. Returns a dict with every
    quantity of the chain (CN, N160, alpha, beta, N160f, CRR75, CM, tau_R, rd, tau_eq, FS) and
    `verdict`: "safe", "triggered", "not susceptible", "not required" or "beyond curve", with
    `note` saying why when no factor is computed."""
    return _core.liquefaction_tbdy2018(
        float(depth), float(N), float(fines), float(sigma_v), float(sigma_v_eff), float(SDS),
        float(Mw), float(CR), float(CS), float(CB), float(CE), bool(below_water),
        bool(susceptible), float(sigma_v_eff_at_test), float(required))


def boring(tests, *, water_depth, gamma, gamma_sat, SDS, Mw=7.5, gamma_w=9.81, rod_factor=None,
           CS=1.0, CB=1.0, CE=1.0, susceptible=None):
    """A boring in uniform ground: `tests` are (depth, N, fines) tuples; the vertical stresses
    are the weight of gamma above the water table and gamma_sat below it, the pore pressure
    hydrostatic from the water table. `rod_factor(depth)` gives C_R (default: Table 16B.1 on
    a rod one metre longer than the depth); `susceptible(depth)` says whether the soil there can
    liquefy at all (default: every tested layer can -- 16.6.2 and 16.6.4 are the engineer's)."""
    rows = []
    for depth, N, fines in tests:
        z = float(depth)
        above = min(z, water_depth)
        below = max(0.0, z - water_depth)
        sv = gamma * above + gamma_sat * below
        u = gamma_w * below
        cr = rod_factor(z) if rod_factor else rod_factor_tbdy(z + 1.0)
        r = spt(z, N, fines, sv, sv - u, SDS=SDS, Mw=Mw, CR=cr, CS=CS, CB=CB, CE=CE,
                below_water=z > water_depth,
                susceptible=True if susceptible is None else bool(susceptible(z)))
        r["depth"] = z
        r["sigma_v"] = sv
        r["sigma_v_eff"] = sv - u
        rows.append(r)
    return rows
