#!/usr/bin/env python3
"""Reference gate: every verification test declares its oracle, or is known not to.

Measured at adoption (2026-07-30): 59 of 126 test files mentioned a source at
all, informally, with no section, equation or tolerance. That is a memory aid,
not a reference. This gate makes the discipline enforceable, and the ratchet
below makes it adoptable without rewriting the whole suite in one day.

The declaration block (roadmap section 6.2), anywhere in a ``tests/test_*.cpp``:

    // verify: KV-FND-001
    //   oracle:   published_benchmark
    //   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems)
    //   locator:  Section 2.1, smooth rigid strip footing on elastic soil
    //   quantity: footing force F at prescribed settlement 10 mm [kN/m]
    //   expected: 15.15 (Giroud analytic); PLAXIS publishes 15.24
    //   band:     2% vs analytic -- coarse-mesh bias measured at +1.4%

Field rules:
  - ``oracle``   one of: closed_form, independent_path, external_code,
                 published_benchmark.
  - ``source``   precise enough to fetch. Never just an author's surname.
  - ``locator``  a section, equation or table number -- or the closed form
                 stated in full, which is just as checkable. Never only a title.
  - ``quantity`` what is compared, with units.
  - ``expected`` the reference value.
  - ``band``     the accepted deviation *and the reason it is that wide*. The
                 band must match what the test actually asserts.

A test that genuinely verifies nothing external (pure regression or unit
mechanics) declares that instead, with a reason:

    // verify: none -- JSON round-trip regression; no external oracle exists

Everything else is the ratchet: a test file with neither a declaration nor a
``verify: none`` marker must appear in UNDECLARED_AT_ADOPTION below. That list
is frozen -- entries leave when a file gains its declaration, and nothing is
ever added. A new test without a declaration fails the gate on arrival.

The verification matrix and references.bib are GENERATED from these blocks
(``--write``), and the gate fails when the committed files differ from
regeneration -- a hand-edited generated table is treated as drift, because a
table maintained beside the tests is exactly what decayed to 59/126.

Usage:
    python scripts/check_references.py            # gate: validate + ratchet + drift
    python scripts/check_references.py --write    # regenerate matrix + bib
    python scripts/check_references.py --list     # print parsed declarations
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = REPO_ROOT / "tests"
MATRIX_PATH = REPO_ROOT / "docs" / "validation" / "verification-matrix.md"
BIB_PATH = REPO_ROOT / "docs" / "validation" / "references.bib"

ORACLES = ("closed_form", "independent_path", "external_code", "published_benchmark")
CLASSES = ("NUM", "CST", "FND", "EXC", "SLP", "CON", "FLW", "STR", "DYN", "SSI")

CASE_RE = re.compile(r"^KV-(" + "|".join(CLASSES) + r")-\d{3}$")
# Only 'KV-...' or 'none ...' payloads open a block: an English sentence that
# happens to wrap onto a line starting with 'verify:' is prose, not a
# declaration -- but a malformed KV id is still an error, never ignored.
START_RE = re.compile(r"^\s*//\s*verify:\s*((?:KV-|none\b)\S*.*?)\s*$")
FIELD_RE = re.compile(r"^\s*//\s+(oracle|source|locator|quantity|expected|band):\s*(.*\S)\s*$")
CONT_RE = re.compile(r"^\s*//\s+(\S.*?)\s*$")
NONE_RE = re.compile(r"^none\s*(?:--|—)\s*(\S.*)$")

# Test files that predate the gate and do not yet carry a declaration.
# Frozen at adoption, 2026-07-30. Entries LEAVE this list when the file gains
# its block; nothing is ever added. A new test without a declaration fails.
UNDECLARED_AT_ADOPTION: frozenset[str] = frozenset({
    "test_anchor.cpp",
    "test_anchor_mesh_repro.cpp",
    "test_anchor_plastic.cpp",
    "test_assembly.cpp",
    "test_axisym_collapse.cpp",
    "test_axisym_cylinder.cpp",
    "test_axisym_gui.cpp",
    "test_axisym_mc.cpp",
    "test_axisymmetric.cpp",
    "test_bearing_phi.cpp",
    "test_compliant_base.cpp",
    "test_consolidation_gui.cpp",
    "test_consolidation_plastic.cpp",
    "test_constitutive.cpp",
    "test_coupled_flow.cpp",
    "test_coupled_flow_plastic.cpp",
    "test_delaunay.cpp",
    "test_design_code.cpp",
    "test_design_gui.cpp",
    "test_distributed_load.cpp",
    "test_dynamic_gui.cpp",
    "test_dynamics.cpp",
    "test_dynamics_nonlinear.cpp",
    "test_earth_pressure.cpp",
    "test_effective_stress.cpp",
    "test_embedded_beam.cpp",
    "test_embedded_wall.cpp",
    "test_excavation.cpp",
    "test_excavation_wall.cpp",
    "test_general_interface.cpp",
    "test_geogrid.cpp",
    "test_geometry2d.cpp",
    "test_gui_matrix.cpp",
    "test_gui_solve.cpp",
    "test_hardening_soil.cpp",
    "test_hs_axisym.cpp",
    "test_hs_berlin.cpp",
    "test_hs_calibration.cpp",
    "test_hs_cap.cpp",
    "test_hs_fe.cpp",
    "test_hs_fe_oedometer.cpp",
    "test_hs_footing.cpp",
    "test_hs_integrate.cpp",
    "test_hs_shear.cpp",
    "test_hs_two_surface.cpp",
    "test_hs_undrained_bvp.cpp",
    "test_hs_undrained_clay.cpp",
    "test_hs_undrained_triaxial.cpp",
    "test_hssmall.cpp",
    "test_initial_stress_k0.cpp",
    "test_input_audit.cpp",
    "test_interface.cpp",
    "test_interface_gui.cpp",
    "test_k0_layered.cpp",
    "test_k0_slope.cpp",
    "test_material_profile.cpp",
    "test_math_smoke.cpp",
    "test_mesh.cpp",
    "test_mesh_domain.cpp",
    "test_mesh_refine.cpp",
    "test_mesher_pipeline.cpp",
    "test_model_mesh.cpp",
    "test_mohr_coulomb.cpp",
    "test_nonlinear.cpp",
    "test_pardiso_smoke.cpp",
    "test_pardiso_solver.cpp",
    "test_phreatic_gravity.cpp",
    "test_plate.cpp",
    "test_plate_mesh.cpp",
    "test_plate_plastic.cpp",
    "test_plate_soil.cpp",
    "test_plate5.cpp",
    "test_point_location.cpp",
    "test_predicates.cpp",
    "test_project_io.cpp",
    "test_results_io.cpp",
    "test_robustness.cpp",
    "test_safety_gui.cpp",
    "test_seepage.cpp",
    "test_seepage_gui.cpp",
    "test_site_response_benchmark.cpp",
    "test_slope.cpp",
    "test_soft_soil.cpp",
    "test_soft_soil_creep.cpp",
    "test_soft_soil_gui.cpp",
    "test_solve.cpp",
    "test_sparse_matrix.cpp",
    "test_ssi_dynamics.cpp",
    "test_staged_fill.cpp",
    "test_staged_gui.cpp",
    "test_staged_struct_carry.cpp",
    "test_stress.cpp",
    "test_struct_forces.cpp",
    "test_struct_gui.cpp",
    "test_tbdy_seismic.cpp",
    "test_transient_flow.cpp",
    "test_tri15.cpp",
    "test_tri15_pipeline.cpp",
    "test_tri15_prandtl.cpp",
    "test_tri15_traction.cpp",
    "test_tri6.cpp",
    "test_undrained.cpp",
    "test_undrained_global.cpp",
    "test_undrained_gui.cpp",
    "test_unsaturated_flow.cpp",
    "test_wall_k0_excavation.cpp",
    "test_wall_split.cpp",
    "test_water_gui.cpp",
    "test_water_retention.cpp",
})


@dataclass
class Declaration:
    case: str
    file: str
    line: int
    fields: dict[str, str] = field(default_factory=dict)

    @property
    def cls(self) -> str:
        return self.case.split("-")[1]


def parse_file(path: Path) -> tuple[list[Declaration], list[str], list[str]]:
    """Returns (declarations, none-markers' reasons, problems)."""
    decls: list[Declaration] = []
    nones: list[str] = []
    problems: list[str] = []
    rel = path.relative_to(REPO_ROOT).as_posix()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    i = 0
    while i < len(lines):
        m = START_RE.match(lines[i])
        if not m:
            i += 1
            continue
        head = m.group(1)
        none = NONE_RE.match(head)
        if none:
            nones.append(none.group(1))
            i += 1
            continue
        decl = Declaration(case=head, file=rel, line=i + 1)
        if not CASE_RE.match(head):
            problems.append(f"{rel}:{i + 1}: bad case id '{head}' (KV-<CLASS>-<nnn>, "
                            f"classes: {', '.join(CLASSES)}; or 'none -- reason')")
        i += 1
        current: str | None = None
        while i < len(lines):
            fm = FIELD_RE.match(lines[i])
            if fm:
                current = fm.group(1)
                if current in decl.fields:
                    problems.append(f"{rel}:{i + 1}: duplicate field '{current}' in {decl.case}")
                decl.fields[current] = fm.group(2)
                i += 1
                continue
            cm = CONT_RE.match(lines[i])
            if cm and current and not START_RE.match(lines[i]):
                decl.fields[current] += " " + cm.group(1)
                i += 1
                continue
            break
        for key in ("oracle", "source", "locator", "quantity", "expected", "band"):
            if not decl.fields.get(key):
                problems.append(f"{rel}:{decl.line}: {decl.case} lacks '{key}'")
        oracle = decl.fields.get("oracle", "")
        if oracle and oracle not in ORACLES:
            problems.append(f"{rel}:{decl.line}: {decl.case} oracle '{oracle}' not in "
                            f"{{{', '.join(ORACLES)}}}")
        decls.append(decl)
    return decls, nones, problems


def collect() -> tuple[list[Declaration], dict[str, bool], list[str]]:
    """All declarations, per-file declared/none status, and problems."""
    all_decls: list[Declaration] = []
    status: dict[str, bool] = {}
    problems: list[str] = []
    for path in sorted(TESTS_DIR.glob("test_*.cpp")):
        decls, nones, probs = parse_file(path)
        problems.extend(probs)
        all_decls.extend(decls)
        status[path.name] = bool(decls) or bool(nones)
    seen: dict[str, Declaration] = {}
    for d in all_decls:
        if CASE_RE.match(d.case):
            if d.case in seen:
                problems.append(f"{d.file}:{d.line}: duplicate case id {d.case} "
                                f"(also in {seen[d.case].file}:{seen[d.case].line})")
            else:
                seen[d.case] = d
    return all_decls, status, problems


def render_matrix(decls: list[Declaration], undeclared: list[str]) -> str:
    rows = sorted((d for d in decls if CASE_RE.match(d.case)), key=lambda d: d.case)
    out: list[str] = []
    out.append("# Verification matrix")
    out.append("")
    out.append("<!-- GENERATED by scripts/check_references.py --write. Do not edit: the gate")
    out.append("     test_reference_registry fails when this file differs from regeneration.")
    out.append("     The source of truth is the declaration block inside each test. -->")
    out.append("")
    out.append("Each row is generated from a machine-readable declaration (roadmap section 6.2)")
    out.append("inside the test that enforces it, so this table cannot drift from the suite.")
    out.append("The `Input file` column stays empty until the input corpus (section 6.3) lands;")
    out.append("every row is asserted by a ctest regression on every build of every composition")
    out.append("that contains its modules.")
    out.append("")
    out.append("| Case | Class | Quantity | Source and locator | Expected | Band | Oracle | Input file | Test |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for d in rows:
        f = d.fields
        test = Path(d.file).stem
        src = f"{f.get('source', '?')} — {f.get('locator', '?')}"
        out.append(f"| {d.case} | {d.cls} | {f.get('quantity', '?')} | {src} "
                   f"| {f.get('expected', '?')} | {f.get('band', '?')} "
                   f"| {f.get('oracle', '?')} | — | `{test}` |")
    out.append("")
    out.append(f"**Coverage, computed:** {len(rows)} declared case(s) across "
               f"{len({d.file for d in rows})} test file(s); "
               f"{len(undeclared)} test file(s) predate the gate and remain in the frozen "
               f"adoption list inside `scripts/check_references.py` — that number may only fall.")
    out.append("")
    return "\n".join(out) + "\n"


def bib_key(source: str, n: int) -> str:
    words = re.findall(r"[A-Za-z]+", source)
    year = re.search(r"(19|20)\d{2}", source)
    head = (words[0].lower() if words else "ref")
    return f"{head}{year.group(0) if year else ''}{'' if n == 0 else chr(ord('a') + n)}"


def render_bib(decls: list[Declaration]) -> str:
    by_source: dict[str, list[str]] = defaultdict(list)
    for d in sorted((x for x in decls if CASE_RE.match(x.case)), key=lambda x: x.case):
        by_source[d.fields.get("source", "?")].append(d.case)
    out = [
        "% GENERATED by scripts/check_references.py --write. Do not edit.",
        "% One entry per distinct declared source; the note lists the cases citing it.",
        "",
    ]
    used: set[str] = set()
    for source in sorted(by_source):
        n = 0
        key = bib_key(source, n)
        while key in used:
            n += 1
            key = bib_key(source, n)
        used.add(key)
        cases = ", ".join(by_source[source])
        out.append(f"@misc{{{key},")
        out.append(f"  title = {{{source}}},")
        out.append(f"  note  = {{cited by {cases}}},")
        out.append("}")
        out.append("")
    return "\n".join(out) + ("\n" if out[-1] else "")


def main() -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")

    ap = argparse.ArgumentParser(description="Verification-reference gate (roadmap section 6.2).")
    ap.add_argument("--write", action="store_true",
                    help="regenerate docs/validation/verification-matrix.md and references.bib")
    ap.add_argument("--list", action="store_true", help="print parsed declarations and exit")
    args = ap.parse_args()

    decls, status, problems = collect()
    undeclared_now = sorted(name for name, ok in status.items() if not ok)

    if args.list:
        for d in sorted((x for x in decls if CASE_RE.match(x.case)), key=lambda x: x.case):
            print(f"{d.case}  {d.file}:{d.line}  [{d.fields.get('oracle', '?')}] "
                  f"{d.fields.get('quantity', '?')}")
        print(f"\n{len(decls)} declaration(s); {len(undeclared_now)} file(s) undeclared.")
        return 0

    matrix = render_matrix(decls, undeclared_now)
    bib = render_bib(decls)

    if args.write:
        MATRIX_PATH.write_text(matrix, encoding="utf-8", newline="\n")
        BIB_PATH.write_text(bib, encoding="utf-8", newline="\n")
        print(f"check_references: wrote {MATRIX_PATH.name} "
              f"({sum(1 for d in decls if CASE_RE.match(d.case))} row(s)) and {BIB_PATH.name}.")
        return 0

    # Gate mode. 1: block validity.
    failed = list(problems)

    # 2: ratchet -- an undeclared test is either historic (frozen list) or an error.
    for name in undeclared_now:
        if name not in UNDECLARED_AT_ADOPTION:
            failed.append(f"tests/{name}: no declaration block and no 'verify: none -- reason' "
                          f"marker. New tests must declare their oracle (roadmap section 6.2).")
    for name in sorted(UNDECLARED_AT_ADOPTION):
        if name not in status:
            failed.append(f"scripts/check_references.py: '{name}' is in the frozen adoption list "
                          f"but no such test exists -- remove the stale entry.")
        elif status[name]:
            failed.append(f"tests/{name} now declares -- remove it from UNDECLARED_AT_ADOPTION "
                          f"so the ratchet advances.")

    # 3: generated artifacts must equal regeneration.
    for path, want in ((MATRIX_PATH, matrix), (BIB_PATH, bib)):
        have = path.read_text(encoding="utf-8") if path.exists() else None
        if have != want:
            failed.append(f"{path.relative_to(REPO_ROOT).as_posix()} differs from regeneration -- "
                          f"run: python scripts/check_references.py --write")

    if failed:
        print(f"check_references: {len(failed)} problem(s)")
        for item in failed:
            print("  " + item)
        return 1

    n_cases = sum(1 for d in decls if CASE_RE.match(d.case))
    print(f"check_references: OK - {n_cases} declared case(s), "
          f"{len(undeclared_now)} file(s) still in the frozen adoption list.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
