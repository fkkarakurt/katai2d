"""The schema-coverage gate: every .k2d key is reachable from Python.

check_k2d_spec.py already locks the format's four descriptions together
(writer = reader = spec = JSON Schema). This test welds the FIFTH surface,
the Python bindings, onto the same set: every key the writer emits must be

  * REACHABLE  -- the corresponding bound object has the attribute,
  * FAITHFUL   -- reading it gives exactly the value the file carries,
  * WRITABLE   -- assigning it back is accepted (no silently read-only field).

The denominator is not hand-maintained: it is extracted from project_io.hpp
by importing check_k2d_spec's own writer_reader_keys(), so a field added to
the schema lands in this test's universe the moment the writer emits it --
and fails here until the binding (and, if the name differs, the KEY_TO_ATTR
row) exists. That is the drift this gate makes impossible to commit silently.

Usage: schema_coverage.py <repo-root>
"""

import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(sys.argv[1]).resolve()

import katai
from katai import _core as core

# ---------------------------------------------------------------- denominator --
# Import the spec gate as a module and reuse its extraction verbatim: one regex,
# one source of truth. (Its main() is __main__-guarded; importing runs nothing.)
_spec = importlib.util.spec_from_file_location(
    "check_k2d_spec", REPO / "scripts" / "check_k2d_spec.py")
_gate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gate)
WRITTEN, _ = _gate.writer_reader_keys(
    (REPO / "kernel/io/include/katai/io/project_io.hpp").read_text(encoding="utf-8-sig"))
assert len(WRITTEN) > 100, f"writer key extraction looks broken ({len(WRITTEN)} keys)"

# File keys and Python attribute names coincide except where the file format
# abbreviates (the file is the older contract; the bindings follow the C++
# member names, which spell things out). A NEW schema field whose file key
# differs from its member name must add its row here -- the failure message
# below says so.
KEY_TO_ATTR = {
    "steps": "time_steps",
    "design": "design_approach",
    "seiswave": "seismic_wave",
    "seisamp": "seismic_amp",
    "seisfreq": "seismic_freq",
    "damp": "damping_ratio",
    "rayf1": "rayleigh_f1",
    "rayf2": "rayleigh_f2",
    "seisff": "seismic_free_field",
    "dynnl": "dynamic_nonlinear",
    "seiscb": "seismic_compliant_base",
    "ec8on": "ec8_enabled",
    "ec8agr": "ec8_agr",
    "ec8gi": "ec8_gamma",
    "ec8gnd": "ec8_ground",
    "ec8typ": "ec8_type",
    "rec": "accel_record",
    "recdt": "record_dt",
    "tbdyss": "tbdy_ss",
    "tbdys1": "tbdy_s1",
    "siteclass": "site_class",
    "mstage": "sum_mstage",
    "ignoreund": "ignore_undrained",
    "tol": "tolerance",
    "loadsteps": "load_steps",
    "maxiter": "max_iterations",
    "lamstar": "lam_star",
    "kapstar": "kap_star",
    "mustar": "mu_star",
    "poly": "poly_active",
    "struct": "struct_active",
    "load": "load_active",
    "disp": "disp_active",
}

# ---------------------------------------------------------- the coverage model --
# One of everything, so the writer emits every key it knows. Values stay at
# their C++ defaults wherever possible: the point is reachability, not physics,
# and defaults are exactly what a forgotten binding would silently lose.


def coverage_project():
    pr = core.Project()
    pr.name = "schema coverage"
    pr.has_water = True
    pr.wx = [0.0, 10.0]
    pr.wy = [4.0, 4.0]

    soil = core.Material()
    soil.name = "Soil"
    pr.materials = [soil]
    plate = core.PlateMaterial()
    plate.name = "Wall"
    pr.plates = [plate]
    anchor = core.AnchorMaterial()
    anchor.name = "Anchor"
    pr.anchors = [anchor]
    geogrid = core.GeogridMaterial()
    geogrid.name = "Geogrid"
    pr.geogrids = [geogrid]
    beam = core.EmbeddedBeamMaterial()
    beam.name = "Pile"
    pr.embedded = [beam]

    poly = core.SoilPolygon()
    poly.name = "Region"
    poly.material = 0
    poly.x = [0.0, 10.0, 10.0, 0.0]
    poly.y = [0.0, 0.0, 5.0, 5.0]
    poly.edge_bc = [core.BCType(0), core.BCType(1), core.BCType(2), core.BCType(1)]
    poly.edge_flow = [core.FlowBCType(0), core.FlowBCType(1), core.FlowBCType(2), core.FlowBCType(1)]
    poly.edge_head = [0.0, 0.0, 4.0, 0.0]
    pr.polygons = [poly]

    wall = core.StructElement()
    wall.kind = core.StructKind(0)
    wall.name = "Wall"
    wall.x1, wall.y1, wall.x2, wall.y2 = 5.0, 5.0, 5.0, 1.0
    pr.structs = [wall]

    strip = core.Load()
    strip.kind = core.LoadKind(0)
    strip.name = "Strip"
    strip.x1, strip.y1, strip.x2, strip.y2 = 2.0, 5.0, 4.0, 5.0
    pr.loads = [strip]

    footing = core.PrescribedDisp()
    footing.name = "Footing"
    footing.x1, footing.y1, footing.x2, footing.y2 = 4.0, 5.0, 6.0, 5.0
    footing.set_uy = True
    footing.uy = -0.01
    pr.disps = [footing]

    # The writer emits "rec"/"recdt" only when a record exists -- give the
    # staged phase one, so the coverage universe includes them.
    ph = core.Phase()
    ph.name = "Shake"
    ph.accel_record = [0.0, 0.5, -0.5]
    ph.record_dt = 0.02
    # Same rule for the numerical controls: written only when set, so the coverage
    # phase sets them.
    ph.sum_mstage = 0.5
    ph.ignore_undrained = True
    ph.tolerance = 1e-5
    ph.load_steps = 12
    ph.max_iterations = 40
    ph.poly_active = [1]
    ph.struct_active = [1]
    ph.load_active = [1]
    ph.disp_active = [1]
    pr.phases = [ph]
    return pr


# ------------------------------------------------------------------- the walk --
COVERED = set()
PROBLEMS = []


def as_number(v):
    """One comparable number for float/int/bool and arithmetic enums."""
    try:
        return float(v)
    except TypeError:
        return float(int(v))


def check_leaf(jval, pyval, path):
    if isinstance(jval, str):
        if pyval != jval:
            PROBLEMS.append(f"{path}: file says {jval!r}, Python reads {pyval!r}")
    elif isinstance(jval, list):
        got = list(pyval)
        if len(got) != len(jval) or any(as_number(a) != as_number(b) for a, b in zip(jval, got)):
            PROBLEMS.append(f"{path}: file says {jval}, Python reads {got}")
    elif as_number(pyval) != as_number(jval):
        PROBLEMS.append(f"{path}: file says {jval}, Python reads {pyval}")


def walk(jobj, pyobj, path):
    for key, jval in jobj.items():
        COVERED.add(key)
        if key == "katai2d":   # file envelope, not a Project member
            if jval != core.PROJECT_FILE_VERSION:
                PROBLEMS.append(f"{path}.katai2d: {jval} != PROJECT_FILE_VERSION")
            continue
        attr = KEY_TO_ATTR.get(key, key)
        if not hasattr(pyobj, attr):
            PROBLEMS.append(
                f"{path}.{key}: no Python attribute {type(pyobj).__name__}.{attr} "
                f"(bind it; if the file key abbreviates, add a KEY_TO_ATTR row)")
            continue
        pyval = getattr(pyobj, attr)
        if isinstance(jval, dict):
            walk(jval, pyval, f"{path}.{key}")
        elif isinstance(jval, list) and jval and isinstance(jval[0], dict):
            if len(pyval) != len(jval):
                PROBLEMS.append(f"{path}.{key}: {len(jval)} in file, {len(pyval)} in Python")
                continue
            for i, (jv, pv) in enumerate(zip(jval, pyval)):
                walk(jv, pv, f"{path}.{key}[{i}]")
        else:
            check_leaf(jval, pyval, f"{path}.{key}")
            # Writable: assigning the value a LEAF already has must be accepted.
            # Only leaves are probed -- writing a container attribute back while
            # the walk still holds its children reassigns the underlying vector
            # and leaves those children dangling (measured: garbage reads plus
            # heap corruption 0xc0000374 at exit). Container writability is
            # already proven by coverage_project(), which assigns every one.
            try:
                setattr(pyobj, attr, pyval)
            except Exception as e:  # noqa: BLE001 -- any refusal is the finding
                PROBLEMS.append(f"{path}.{key}: {type(pyobj).__name__}.{attr} is not writable ({e})")


pr = coverage_project()
walk(json.loads(core.project_to_json(pr)), pr, "$")

missing = WRITTEN - COVERED
for key in sorted(missing):
    PROBLEMS.append(
        f"writer key `{key}` never appeared in the coverage project's JSON -- "
        f"extend coverage_project() so the writer emits it")
stray = COVERED - WRITTEN
for key in sorted(stray):
    PROBLEMS.append(f"walked key `{key}` is not a writer key (extraction drift?)")

if PROBLEMS:
    for p in PROBLEMS:
        print(f"  {p}")
    print(f"schema_coverage: {len(PROBLEMS)} problem(s)")
    sys.exit(1)

print(f"ok: all {len(WRITTEN)} .k2d keys are reachable, faithful and writable from Python")
