"""Slope stability with KATAI 2D, end to end -- a first script.

The problem: a homogeneous 1:2 slope on a foundation layer, the simple slope of
the referee problems in Giam & Donald (1989), "Example problems for testing soil
slope stability programs", Monash University report 8/1989. Referee factor of
safety: 1.00. Strength reduction by finite elements after Griffiths & Lane (1999).

Units everywhere: kN, m, day. Run with the built package on PYTHONPATH::

    python slope_stability.py
"""
import katai

# 1. A project: plane strain, quadratic elements, ~3 m target mesh size.
prj = katai.Project("Griffiths & Lane slope", mesh_size=3.0, auto_refine=False)

# 2. One soil, Mohr-Coulomb: E [kN/m2], c [kN/m2], phi [deg], gamma [kN/m3].
soil = prj.materials.mohr_coulomb("Clayey sand", E=1.0e5, nu=0.3,
                                  c=3.0, phi=19.6, gamma=20.2)

# 3. The geometry, drawn counter-clockwise; one fixity name per edge
#    (base fixed, far sides on rollers, ground surface free).
prj.geometry.polygon(
    [(20, 20), (70, 20), (70, 35), (50, 35), (30, 25), (20, 25)],
    material=soil,
    fix=["full", "horizontal", "free", "free", "free", "horizontal"])

# 4. The analysis: a phi-c reduction of the gravity state = the slope's FoS.
prj.initial(procedure="safety")

# 5. Run -- the same job the GUI's Calculate button and `katai solve` submit.
job = prj.run(on_phase=lambda i, n, name: print(f"  phase {i + 1}/{n}: {name}"))
res = job.results()[-1]

sign = ">" if res.fos_is_lower_bound else "="
print(f"Factor of safety {sign} {res.fos:.3f}   (referee value: 1.00)")
print(f"Mechanism max displacement: {res.max_disp:.3e} m (unit-less shape)")

# The same analysis, saved for the GUI or the command line:
#   prj.save("slope.k2d")        ->  katai solve slope.k2d

# The corpus band for this exact case, on this exact mesh (KV-SLP-001).
if abs(res.fos - 1.00) > 0.08 * 1.00:
    raise SystemExit("FoS drifted outside the corpus band -- investigate")
