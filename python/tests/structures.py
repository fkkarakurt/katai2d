"""Structural elements from the script, and the output you get back.

Until now the pleasant surface stopped exactly where real engineering starts: a
script could build soil, water, loads and phases, but not a wall, an anchor, a
geogrid or a pile -- for those you had to drop to the raw bindings. This test
covers the layer that closes that, and the summary layer that makes the answer
readable without hand-rolling a numpy scan.

Why this case is checked by SOLVING rather than by byte-identity, unlike
dsl_corpus.py: KV-STR-005's C++ builder writes the "everything is active"
activation vectors explicitly, and the DSL writes a vector only for a class a
phase actually touches. An omitted vector MEANS all-active (Phase::active_struct),
so the two files describe the same run and differ only in redundancy -- and making
the DSL emit them would break the twelve cases that ARE byte-identical today. So
the property asserted here is the one that matters: the same answer, bit for bit.
"""
import sys

import katai

CORPUS = sys.argv[1]
failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


def geogrid_case(*, elastoplastic=False, Np=3.0):
    """KV-STR-005 as a script: a stretched elastic block with a full-width sheet."""
    p = katai.Project("KV-STR-005 geogrid axial force and tension cut-off", mesh_size=0.5)
    fill = p.materials.linear_elastic("Fill", E=10000, nu=0.3, c=1, phi=30,
                                      gamma=0, gamma_sat=0)
    p.geometry.polygon([(0, 0), (10, 0), (10, 5), (0, 5)], material=fill, name="Fill",
                       fix=["vertical", "free", "free", "horizontal"])
    p.structures.geogrid((0, 2), (10, 2), EA=5000, Np=Np, elastoplastic=elastoplastic,
                         name="Reinforcement")
    d = p.displacements.line((10, 0), (10, 5), ux=0.004, name="Stretch")
    p.initial(exclude=[d])
    p.phases.plastic("Stretch", activate=[d])
    return p


# ---------------------------------------------------------------- the same run --
script = geogrid_case()
from_file = katai.load_project(f"{CORPUS}/kv-str-005-geogrid-tension.k2d")[0]

a = katai.run(script.build()).results()[-1]
b = katai.run(from_file).results()[-1]
check(a.ok and b.ok, "both the scripted and the checked-in geogrid case solve")
check(a.max_disp == b.max_disp,
      f"the script and the file give the SAME answer, bit for bit ({a.max_disp:.12e} m)")
check(len(a.struct_forces) == 1 and len(b.struct_forces) == 1,
      "one structural element reports its forces on both paths")
if a.struct_forces and b.struct_forces:
    check(a.struct_forces[0].max_N == b.struct_forces[0].max_N,
          f"and the same geogrid force ({a.struct_forces[0].max_N:.9f} kN/m)")
    # The case's own oracle: EA * eps = 5000 * 0.004/10 = 2.0 kN/m.
    check(abs(a.struct_forces[0].max_N - 2.0) < 2e-3,
          "which is EA*eps = 2.000 kN/m, the closed form KV-STR-005 verifies")

# A capacity the engine never reads is worse than no capacity: the driver applies a
# geogrid cap only when the element is marked elastoplastic. Giving Np must therefore
# ARM it, or `Np=3` would be a silent no-op in a script that believed it set a limit.
armed = geogrid_case(elastoplastic=None, Np=1.0).build()
check(armed.geogrids[0].elastoplastic and armed.geogrids[0].Np == 1.0,
      "a stated tension cap arms itself (Np= is not a silent no-op)")
dormant = geogrid_case(elastoplastic=False, Np=1.0).build()
check(not dormant.geogrids[0].elastoplastic,
      "...and elastoplastic=False can still state a dormant capacity, as the corpus file does")
capped = katai.run(geogrid_case(elastoplastic=True, Np=1.0).build()).results()[-1]
check(capped.struct_forces[0].max_N <= 1.0 + 1e-9,
      f"an armed cap actually bites: N = {capped.struct_forces[0].max_N:.6f} <= 1.0 kN/m")

# ------------------------------------------------------- every kind round-trips --
p = katai.Project("All five kinds", mesh_size=2.0, auto_refine=False)
soil = p.materials.mohr_coulomb("Sand", E=3e4, nu=0.3, c=5, phi=32, gamma=19)
p.geometry.rectangle(0, 0, 30, 20, material=soil, name="Ground")
p.structures.plate((12, 20), (12, 6), EA=1.2e7, EI=1.0e5, w=5.0, nu=0.15,
                   interfaces="both", barrier="impermeable", name="Wall")
p.structures.anchor((12, 18), (20, 14), EA=2.0e5, spacing=2.5, prestress=300,
                    Fmax_tens=800, name="Anchor")
p.structures.geogrid((0, 8), (10, 8), EA=5000, Np=3, name="Sheet")
p.structures.pile((22, 20), (22, 8), E=3e7, diameter=0.6, spacing=2.5,
                  Tskin_max=50, Fmax_base=200, connection="free", name="Pile")
p.structures.interface((5, 20), (5, 10), name="Slip line")
built = p.build()
check(len(built.structs) == 5, "all five structural kinds reach the schema")
check([int(s.kind) for s in built.structs] == [0, 1, 2, 3, 4],
      "...each as its own kind, in the order the script drew them")
check(built.plates[0].EA == 1.2e7 and built.plates[0].nu == 0.15,
      "plate properties carry through")
check(built.structs[0].iface_pos and built.structs[0].iface_neg,
      "interfaces='both' puts one on each side of the wall")
check(built.structs[0].flow_barrier == 1, "a wall can be made a groundwater screen")
check(built.anchors[0].prestress == 300 and built.anchors[0].Lspacing == 2.5,
      "anchor lock-off force and out-of-plane spacing carry through")
check(built.anchors[0].elastoplastic and built.anchors[0].Fmax_tens == 800,
      "a stated anchor capacity arms itself too")
check(built.embedded[0].diameter == 0.6 and built.embedded[0].Tskin_max == 50,
      "pile geometry and shaft capacity carry through")
check(built.structs[3].conn == 1, "connection='free' reaches the schema (hinged is the default)")
check(built.structs[4].material == -1, "a bare interface takes the adjacent soil, not a material")
check(katai.validate_project(built).ok(), "the five-element project satisfies the input contract")

# The file is the contract: everything above must survive a save/load round trip.
text = katai.project_to_json(built)
back, _ = katai.project_from_json(text)
check(katai.project_to_json(back) == text, "structures round-trip through the .k2d text exactly")

# The extent has to contain the elements, or the mesher never sees them: the wall here
# toes at y = 6 inside the block, but a pile driven BELOW the soil must still fit.
deep = katai.Project("Deep pile", mesh_size=2.0)
s2 = deep.materials.mohr_coulomb("Sand", E=3e4, nu=0.3, c=5, phi=32, gamma=19)
deep.geometry.rectangle(0, 0, 20, 10, material=s2)
deep.structures.pile((10, 10), (10, -4), E=3e7, diameter=0.5, spacing=2.0)
check(deep.build().y_min == -4, "the model extent grows to contain a structural element")

# ---------------------------------------------------------------- phase toggles --
p2 = katai.Project("Staged", mesh_size=2.0, auto_refine=False)
s3 = p2.materials.mohr_coulomb("Sand", E=3e4, nu=0.3, c=5, phi=32, gamma=19)
p2.geometry.rectangle(0, 0, 20, 10, material=s3)
strut = p2.structures.anchor((0, 8), (20, 8), EA=1e5, spacing=2.0, name="Strut")
p2.initial(exclude=[strut])
p2.phases.plastic("Install strut", activate=[strut])
b2 = p2.build()
check(list(b2.initial.struct_active) == [0], "a structure can be absent from the initial phase")
check(list(b2.phases[0].struct_active) == [1], "...and installed by a later phase, by handle")

# ------------------------------------------------------------- readable output --
job = katai.run(geogrid_case().build())
text = katai.summary(job, script)
check("Stretch" in text, "the summary names the phases the script named")
check("Initial phase" in text, "...including the initial one")
check("Reinforcement" in text and "geogrid" in text,
      "the summary reports the structural element by name and kind")
check("extremes" in text and "u_x" in text, "the summary tabulates the nodal extremes")
check("converged" in text, "the summary states the convergence outcome")
ex = katai.extremes(job.results()[-1])
check(set(("|u|", "u_x", "u_y", "p_w")) <= set(ex), "extremes() covers the usual fields")
lo, plo, hi, phi = ex["u_x"]
check(hi > lo and plo is not None and len(phi) == 2,
      "an extreme comes with the place it occurs, which is half the finding")
flat_lo, flat_at, flat_hi, _ = ex["p_w"]
check(flat_lo == flat_hi and flat_at is None,
      "a field that is the same everywhere reports no location rather than inventing one")

print(f"\n{failures} check(s) failed" if failures else "\nall checks passed")
sys.exit(1 if failures else 0)
