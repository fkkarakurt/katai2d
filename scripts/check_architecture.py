#!/usr/bin/env python3
"""Architecture gate: fail the build when the layer contract is violated.

The module layout only stays true if something enforces it. This gate reads the
include graph out of the source tree -- no compilation, so it costs a second
rather than an hour -- and refuses three specific classes of decay:

  dag        Module dependencies must form a directed acyclic graph and may only
             point downward. ``mesh`` may include ``geometry``; ``geometry`` may
             not include ``mesh``. The application may include the engine; the
             engine may never include the application.

  api        Only ``include/katai/**`` is public. One module reaching into
             another module's ``src/`` or into a private header is a violation
             even when it compiles.

  physics    No constitutive, assembly or analysis code above the engine. This
             is the rule that keeps the split real, so it is checked by symbol
             definition rather than by file name.

  config     Optional capabilities are separate targets, not ``#ifdef``s inside
             shared code. Conditional compilation scattered through the tree
             multiplies the number of programs that exist while testing only one
             of them, which is how a configuration nobody built becomes a
             configuration nobody can trust. ``KATAI_WITH_*`` is therefore
             allowed only at a composition root.

Each check reports every violation it finds and the gate exits non-zero if any
check fails. Known, accepted violations live in ``ARCHITECTURE_EXCEPTIONS`` with
the reason and the stage that removes them: an exception is a debt with an
owner, not a permission.

It also answers a structural question that comes up whenever a module might be
released on its own: what would have to go with it. A module cannot be published
without its dependency closure, and it should not be published without the tests
that cover it, because an unverifiable module is worth nothing to a reader. Both
are computed from the tree rather than remembered.

Usage:
    python scripts/check_architecture.py            # all checks, exit 1 on any
    python scripts/check_architecture.py --check dag
    python scripts/check_architecture.py --list     # print the resolved graph
    python scripts/check_architecture.py --closure  # what each module drags with it
    python scripts/check_architecture.py --closure materials
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

SUFFIXES = {".cpp", ".hpp", ".h", ".cc", ".cxx", ".inl", ".ipp"}

# Layer number per module. A module may depend on its own layer or below.
# 1 engine · 2 execution · 3 published API · 4 front ends.
LAYERS: dict[str, int] = {
    "math": 1,
    "linsolve": 1,
    "geometry": 1,
    "mesh": 1,
    "materials": 1,
    "fem": 1,
    "analysis": 1,
    "model": 1,
    "io": 1,
    "jobs": 2,
    "api": 3,
    "render": 4,
    "app": 4,
    "cli": 4,
    "python": 4,
}

# Where each module's sources live, relative to the repository root.
MODULE_ROOTS: dict[str, str] = {
    "math": "kernel/math",
    "linsolve": "kernel/linsolve",
    "geometry": "kernel/geometry",
    "mesh": "kernel/mesh",
    "materials": "kernel/materials",
    "fem": "kernel/fem",
    "analysis": "kernel/analysis",
    "model": "kernel/model",
    "io": "kernel/io",
    "jobs": "kernel/jobs",
    "api": "kernel/api",
    "render": "studio/render",
    "app": "studio/app",
    "cli": "cli",
    "python": "python",
}

# Modules that may legitimately be depended on by anything (leaf utilities).
# Empty on purpose: every dependency is directional and stated above.
UNIVERSAL: set[str] = set()

INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]katai/([a-z0-9_]+)/([^>"]+)[>"]', re.MULTILINE)
PRIVATE_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]*(?:/src/|(?:^|/)detail/)[^>"]*)[>"]', re.MULTILINE)

# Symbols that constitute physics. A *definition* of one of these above the
# engine means the boundary has leaked, whatever the file is called. Calling them
# is exactly what a front end is supposed to do, so the check below separates a
# definition from a call rather than matching the name alone.
PHYSICS_SIGNATURES = (
    re.compile(r"(?<![:.>\w])(?:integrate_point|return_mapping|consistent_tangent|yield_function)\s*\("),
    re.compile(r"(?<![:.>\w])assemble_(?:stiffness|global|element)\w*\s*\("),
    re.compile(r"(?<![:.>\w])element_stiffness\w*\s*\("),
)

# What may sit between a signature's closing parenthesis and its body.
TRAILING_SPECIFIERS = re.compile(r"\s*(?:const|noexcept|override|final|&|&&)\s*")

# A composition root selects what the program is made of, so it is the one place
# a build-configuration macro belongs.
COMPOSITION_ROOTS = (
    "studio/app/main.cpp",
    "cli/main.cpp",
)

CONFIG_MACRO_RE = re.compile(r"\bKATAI_WITH_[A-Z0-9_]+")

# Accepted violations: (module, path substring, reason, removed by).
# All three application-layer driver exceptions were discharged in Stage D:
# build_problem/build_mesh/build_flow relocated to kernel/jobs (layer 2) as
# driver.hpp / mesh_builder.hpp / flow_driver.hpp. The mesh builder went to jobs,
# not to the mesh module as once planned, because it resolves katai::model data
# and the mesh module deliberately stays schema-free.
ARCHITECTURE_EXCEPTIONS: tuple[tuple[str, str, str, str], ...] = ()


def rel(path: Path) -> str:
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def excepted(module: str, path: Path) -> tuple[str, str] | None:
    where = rel(path)
    for exc_module, needle, reason, removed_by in ARCHITECTURE_EXCEPTIONS:
        if exc_module == module and needle in where:
            return reason, removed_by
    return None


def owning_module(path: Path) -> str | None:
    where = rel(path)
    best: tuple[int, str] | None = None
    for module, root in MODULE_ROOTS.items():
        if where.startswith(root + "/"):
            if best is None or len(root) > best[0]:
                best = (len(root), module)
    return best[1] if best else None


def module_sources() -> dict[str, list[Path]]:
    found: dict[str, list[Path]] = defaultdict(list)
    for module, root in MODULE_ROOTS.items():
        base = REPO_ROOT / root
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix in SUFFIXES:
                if owning_module(path) == module:
                    found[module].append(path)
    return found


def read(path: Path) -> str:
    try:
        with open(path, "r", encoding="utf-8", newline="") as fh:
            return fh.read()
    except (UnicodeDecodeError, OSError):
        return ""


def check_dag(sources: dict[str, list[Path]]) -> list[str]:
    """Dependencies point downward only, and never cycle inside a layer."""
    problems: list[str] = []
    edges: dict[str, set[str]] = defaultdict(set)

    for module, paths in sources.items():
        for path in paths:
            for dep, _ in INCLUDE_RE.findall(read(path)):
                if dep == module or dep in UNIVERSAL:
                    continue
                if dep not in LAYERS:
                    problems.append(f"{rel(path)}: includes unknown module 'katai/{dep}/'")
                    continue
                edges[module].add(dep)
                if LAYERS[dep] > LAYERS[module]:
                    note = excepted(module, path)
                    if note:
                        continue
                    problems.append(
                        f"{rel(path)}: layer {LAYERS[module]} module '{module}' includes "
                        f"layer {LAYERS[dep]} module '{dep}' - dependencies point downward only"
                    )

    # Cycles can only occur between modules on the same layer; find them there.
    colour: dict[str, int] = defaultdict(int)
    stack: list[str] = []

    def visit(node: str) -> None:
        colour[node] = 1
        stack.append(node)
        for nxt in sorted(edges.get(node, ())):
            if LAYERS.get(nxt) != LAYERS.get(node):
                continue
            if colour[nxt] == 1:
                cut = stack[stack.index(nxt):]
                problems.append("include cycle: " + " -> ".join(cut + [nxt]))
            elif colour[nxt] == 0:
                visit(nxt)
        stack.pop()
        colour[node] = 2

    for module in sorted(edges):
        if colour[module] == 0:
            visit(module)
    return problems


def check_api(sources: dict[str, list[Path]]) -> list[str]:
    """Only include/katai/** is public; nobody reaches into another module's src."""
    problems: list[str] = []
    for module, paths in sources.items():
        for path in paths:
            text = read(path)
            for bad in PRIVATE_INCLUDE_RE.findall(text):
                problems.append(f"{rel(path)}: includes private path '{bad}'")
            for dep, header in INCLUDE_RE.findall(text):
                if dep == module:
                    continue
                root = MODULE_ROOTS.get(dep)
                if root and not (REPO_ROOT / root / "include" / "katai" / dep / header).exists():
                    problems.append(
                        f"{rel(path)}: includes <katai/{dep}/{header}>, which is not a public "
                        f"header of module '{dep}'"
                    )
    return problems


def is_definition(text: str, open_paren: int) -> bool:
    """True when the parameter list at ``open_paren`` is followed by a body.

    A definition reads ``name(args) { ... }``; a call reads ``name(args);`` or
    appears inside a larger expression. Walking the parenthesis depth is exact
    where a line-based guess is not, and it costs nothing at this scale.
    """
    depth = 0
    i = open_paren
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                rest = text[i + 1 :]
                while True:
                    trimmed = TRAILING_SPECIFIERS.match(rest)
                    if not trimmed or trimmed.end() == 0:
                        break
                    rest = rest[trimmed.end() :]
                return rest.lstrip().startswith("{")
        elif ch in "\"'":
            quote = ch
            i += 1
            while i < n and text[i] != quote:
                i += 2 if text[i] == "\\" else 1
        i += 1
    return False


def check_physics(sources: dict[str, list[Path]]) -> list[str]:
    """No constitutive, assembly or analysis definitions above the engine."""
    problems: list[str] = []
    for module, paths in sources.items():
        if LAYERS.get(module, 1) == 1:
            continue
        for path in paths:
            if excepted(module, path):
                continue
            text = read(path)
            for signature in PHYSICS_SIGNATURES:
                for match in signature.finditer(text):
                    open_paren = text.index("(", match.start())
                    if not is_definition(text, open_paren):
                        continue
                    line = text[: match.start()].count("\n") + 1
                    problems.append(
                        f"{rel(path)}:{line}: physics defined above the engine "
                        f"('{match.group(0).rstrip('(').strip()}' in layer "
                        f"{LAYERS[module]} module '{module}')"
                    )
                    break
    return problems


def check_config(sources: dict[str, list[Path]]) -> list[str]:
    """Build-configuration macros belong at a composition root and nowhere else."""
    problems: list[str] = []
    for module, paths in sources.items():
        for path in paths:
            where = rel(path)
            if any(where.endswith(root) for root in COMPOSITION_ROOTS):
                continue
            if excepted(module, path):
                continue
            text = read(path)
            match = CONFIG_MACRO_RE.search(text)
            if match:
                line = text[: match.start()].count("\n") + 1
                problems.append(
                    f"{where}:{line}: '{match.group(0)}' outside a composition root - an "
                    f"optional capability is a separate target, not a conditional"
                )
    return problems


# Front ends program against the published facade: layer 3 (api), layer 2 for runs
# and progress (jobs), and their own presentation module (render). A front end that
# includes an engine header directly has stepped around the section 7.3 boundary,
# even though the layer DAG alone would allow the downward edge (Stage D gate:
# "no front end includes an engine internal").
FRONTEND_MODULES = ("app", "cli", "python")
FRONTEND_ALLOWED = {"api", "jobs", "render"}


def check_frontend(sources: dict[str, list[Path]]) -> list[str]:
    """Front ends include only the facade, the jobs layer and render."""
    problems: list[str] = []
    for module in FRONTEND_MODULES:
        for path in sources.get(module, ()):
            text = read(path)
            for dep, header in INCLUDE_RE.findall(text):
                if dep == module or dep in FRONTEND_ALLOWED:
                    continue
                if excepted(module, path):
                    continue
                problems.append(
                    f"{rel(path)}: front end includes engine header 'katai/{dep}/{header}' - "
                    f"front ends program against <katai/api/katai.hpp> (plus jobs and render)"
                )
    return problems


CHECKS = {
    "dag": check_dag,
    "api": check_api,
    "physics": check_physics,
    "config": check_config,
    "frontend": check_frontend,
}

TEST_DECL_RE = re.compile(r"^\s*katai_add_test\(\s*([A-Za-z0-9_]+)((?:\s+katai::[a-z]+)*)", re.MULTILINE)


def module_edges(sources: dict[str, list[Path]]) -> dict[str, set[str]]:
    edges: dict[str, set[str]] = defaultdict(set)
    for module, paths in sources.items():
        for path in paths:
            for dep, _ in INCLUDE_RE.findall(read(path)):
                if dep != module and dep in LAYERS:
                    edges[module].add(dep)
    return edges


def closure_of(module: str, edges: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    stack = [module]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        stack.extend(edges.get(current, ()))
    return seen


def tests_linking(modules: set[str]) -> list[str]:
    """Tests whose declared link targets fall inside a set of modules."""
    path = REPO_ROOT / "tests" / "CMakeLists.txt"
    if not path.exists():
        return []
    found: list[str] = []
    for name, libs in TEST_DECL_RE.findall(read(path)):
        linked = {lib.split("::")[1] for lib in libs.split() if "::" in lib}
        # katai_add_test defaults to katai::math when no target is named.
        if not linked:
            linked = {"math"}
        if linked <= modules:
            found.append(name)
    return found


def report_closure(sources: dict[str, list[Path]], only: str | None) -> int:
    edges = module_edges(sources)
    present = sorted(m for m in sources if sources[m])
    if only and only not in present:
        print(f"unknown or empty module '{only}'; present: {', '.join(present)}")
        return 1

    for module in ([only] if only else present):
        needed = closure_of(module, edges) & set(present)
        files = sum(len(sources[m]) for m in needed)
        upward = sorted(m for m in needed if LAYERS.get(m, 1) > LAYERS.get(module, 1))
        covering = tests_linking(needed)

        print(f"\n{module}  (layer {LAYERS.get(module, '?')})")
        print(f"  must ship with : {', '.join(sorted(needed - {module})) or '(nothing — a leaf)'}")
        print(f"  total           : {files} files in {len(needed)} module(s)")
        print(f"  tests confined to the closure: {len(covering)}"
              + (f" ({', '.join(covering[:6])}{', ...' if len(covering) > 6 else ''})" if covering else ""))
        if upward:
            print(f"  WARNING: closure reaches a higher layer: {', '.join(upward)}")
    print("\nA module is publishable on its own only together with this closure, and it is only "
          "\nverifiable if the tests that cover it travel with it. Test attribution is read from the "
          "\ndeclared link targets, so a test that includes a header-only module without linking it is "
          "\nnot counted -- treat the test column as a lower bound.")
    return 0


def main() -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")

    ap = argparse.ArgumentParser(description="Fail on violations of the KATAI layer contract.")
    ap.add_argument("--check", choices=sorted(CHECKS), action="append",
                    help="run only the named check (repeatable; default: all)")
    ap.add_argument("--list", action="store_true", help="print the resolved module graph and exit")
    ap.add_argument("--closure", nargs="?", const="", metavar="MODULE",
                    help="what a module would have to be published with (all modules if omitted)")
    args = ap.parse_args()

    sources = module_sources()

    if args.closure is not None:
        return report_closure(sources, args.closure or None)

    if args.list:
        for module in sorted(sources, key=lambda m: (LAYERS.get(m, 9), m)):
            deps: set[str] = set()
            for path in sources[module]:
                deps.update(dep for dep, _ in INCLUDE_RE.findall(read(path)) if dep != module)
            print(f"L{LAYERS.get(module, '?')} {module:10s} "
                  f"{len(sources[module]):3d} files -> {', '.join(sorted(deps)) or '(none)'}")
        return 0

    selected = args.check or sorted(CHECKS)
    failed = 0
    for name in selected:
        problems = CHECKS[name](sources)
        if problems:
            failed += 1
            print(f"check_architecture [{name}]: {len(problems)} violation(s)")
            for item in problems:
                print("  " + item)
        else:
            print(f"check_architecture [{name}]: OK")

    if ARCHITECTURE_EXCEPTIONS:
        print(f"\n{len(ARCHITECTURE_EXCEPTIONS)} accepted exception(s), each owned by a stage:")
        for module, needle, reason, removed_by in ARCHITECTURE_EXCEPTIONS:
            print(f"  {needle}: {reason} - removed by {removed_by}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
