#!/usr/bin/env python3
"""check_k2d_spec: the .k2d format's four descriptions may not drift apart.

The format lives in four places, each with a different audience:

  1. the WRITER   (project_to_json, kernel/io/.../project_io.hpp)  -- the code
  2. the READER   (project_from_json, same header)                 -- the code
  3. the SPEC     (docs/k2d-format.md, field-inventory tables)     -- engineers
  4. the SCHEMA   (docs/k2d.schema.json)                           -- machines

A hand-maintained format document drifts the day someone adds a field and
forgets one of the four. This gate makes that impossible to do silently:

  * the KEY SET written by the writer must equal the set the reader reads
    (an asymmetry is a field that saves but never loads, or vice versa);
  * every writer key must appear in the spec's field-inventory tables, and
    the spec may not document a key the writer no longer emits;
  * the JSON Schema must declare exactly the same keys;
  * the format versions stated in the spec must equal kProjectFileVersion /
    kResultsFileVersion in the headers, and the schema's katai2d maximum;
  * the schema's file-stable enum bounds (initial procedure, soil model,
    drainage, phase type, waveform, design approach, element kind, load kind,
    edge BC, flow BC)
    must equal the enum sizes in kernel/model/.../project.hpp -- appending an
    enum value without widening the schema fails here, not in a user's hands.

Pure source analysis, no compilation (like check_architecture.py): parsing two
headers, one markdown file and one JSON file costs milliseconds.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PROJECT_IO = REPO / "kernel/io/include/katai/io/project_io.hpp"
RESULTS_IO = REPO / "kernel/io/include/katai/io/results_io.hpp"
PROJECT_HPP = REPO / "kernel/model/include/katai/model/project.hpp"
SPEC = REPO / "docs/k2d-format.md"
SCHEMA = REPO / "docs/k2d.schema.json"

# The writer emits keys through wfield/wkey/warr (and wphase's key argument);
# wcolor always writes the literal key "color". The first argument is pinned to
# the output buffer `o` on purpose: the phases loop calls wphase with a
# throwaway key into a scratch string and strips it, and that call is the one
# w-call whose buffer is not named o.
WRITE_RE = re.compile(r'\bw(?:field|key|arr|phase)\s*\(\s*o\s*,\s*"([^"]+)"')
# The reader looks keys up through the Json accessors and find().
READ_RE = re.compile(r'\.(?:num|flag|str|nums|ints|chars)\s*\(\s*"([^"]+)"')
FIND_RE = re.compile(r'\bfind\s*\(\s*"([^"]+)"')
# Spec inventory rows: a table line whose first cell is a backtick-quoted key.
SPEC_ROW_RE = re.compile(r"^\|\s*`([A-Za-z0-9_]+)`\s*\|")
ENUM_RE = re.compile(r"enum class (\w+)\s*[^{]*\{([^}]*)\}", re.S)


def writer_reader_keys(text: str) -> tuple[set[str], set[str]]:
    writer = set(WRITE_RE.findall(text))
    if "wcolor(" in text:
        writer.add("color")
    reader = set(READ_RE.findall(text)) | set(FIND_RE.findall(text))
    return writer, reader


def spec_keys(text: str) -> set[str]:
    keys: set[str] = set()
    inside = False
    for line in text.splitlines():
        if "field-inventory:begin" in line:
            inside = True
        elif "field-inventory:end" in line:
            inside = False
        elif inside:
            m = SPEC_ROW_RE.match(line)
            if m and m.group(1) not in ("Key",):
                keys.add(m.group(1))
    return keys


def schema_keys(node) -> set[str]:
    keys: set[str] = set()
    if isinstance(node, dict):
        props = node.get("properties")
        if isinstance(props, dict):
            keys.update(props.keys())
        for v in node.values():
            keys.update(schema_keys(v))
    elif isinstance(node, list):
        for v in node:
            keys.update(schema_keys(v))
    return keys


def enum_sizes(text: str) -> dict[str, int]:
    sizes: dict[str, int] = {}
    for name, body in ENUM_RE.findall(text):
        sizes[name] = len([e for e in body.split(",") if e.strip()])
    return sizes


def main() -> int:
    problems: list[str] = []

    # utf-8-sig: identical to utf-8 for BOM-less files, and tolerant of the BOM
    # Windows editors like to prepend -- a BOM must not crash a build gate.
    io_text = PROJECT_IO.read_text(encoding="utf-8-sig")
    res_text = RESULTS_IO.read_text(encoding="utf-8-sig")
    model_text = PROJECT_HPP.read_text(encoding="utf-8-sig")
    spec_text = SPEC.read_text(encoding="utf-8-sig")
    schema = json.loads(SCHEMA.read_text(encoding="utf-8-sig"))

    written, read = writer_reader_keys(io_text)
    documented = spec_keys(spec_text)
    declared = schema_keys(schema)

    def diff(missing: set[str], where: str, verb: str) -> None:
        for key in sorted(missing):
            problems.append(f"{where}: `{key}` {verb}")

    diff(written - read, "project_io.hpp", "is written but never read back")
    diff(read - written, "project_io.hpp", "is read but never written (dead or misspelled key)")
    diff(written - documented, "docs/k2d-format.md", "is emitted by the writer but not in the field inventory")
    diff(documented - written, "docs/k2d-format.md", "is documented but the writer no longer emits it")
    diff(written - declared, "docs/k2d.schema.json", "is emitted by the writer but absent from the schema")
    diff(declared - written, "docs/k2d.schema.json", "is declared but the writer does not emit it")

    # -- version identities ---------------------------------------------------
    def const_int(text: str, name: str, where: str) -> int | None:
        m = re.search(rf"{name}\s*=\s*(\d+)", text)
        if not m:
            problems.append(f"{where}: cannot find {name}")
            return None
        return int(m.group(1))

    def spec_int(pattern: str, label: str) -> int | None:
        m = re.search(pattern, spec_text)
        if not m:
            problems.append(f"docs/k2d-format.md: cannot find the '{label}' statement")
            return None
        return int(m.group(1))

    k2d_code = const_int(io_text, "kProjectFileVersion", "project_io.hpp")
    res_code = const_int(res_text, "kResultsFileVersion", "results_io.hpp")
    k2d_spec = spec_int(r"Current `\.k2d` version: \*\*(\d+)\*\*", "Current .k2d version")
    res_spec = spec_int(r"Current `\.res` version: \*\*(\d+)\*\*", "Current .res version")
    if None not in (k2d_code, k2d_spec) and k2d_code != k2d_spec:
        problems.append(f"docs/k2d-format.md: .k2d version {k2d_spec} != kProjectFileVersion {k2d_code}")
    if None not in (res_code, res_spec) and res_code != res_spec:
        problems.append(f"docs/k2d-format.md: .res version {res_spec} != kResultsFileVersion {res_code}")
    schema_ver = schema.get("properties", {}).get("katai2d", {}).get("maximum")
    if k2d_code is not None and schema_ver != k2d_code:
        problems.append(f"docs/k2d.schema.json: katai2d maximum {schema_ver} != kProjectFileVersion {k2d_code}")

    # -- file-stable enum bounds ----------------------------------------------
    sizes = enum_sizes(model_text)

    def bound(path: list[str]) -> object:
        node: object = schema
        for step in path:
            if not isinstance(node, dict) or step not in node:
                return None
            node = node[step]
        return node

    expectations = [
        (["properties", "initial_procedure", "maximum"], "InitialProcedure"),
        (["$defs", "material", "properties", "model", "maximum"], "SoilModel"),
        (["$defs", "material", "properties", "drainage", "maximum"], "Drainage"),
        (["$defs", "phase", "properties", "type", "maximum"], "PhaseType"),
        (["$defs", "phase", "properties", "seiswave", "maximum"], "SeismicWave"),
        (["$defs", "phase", "properties", "design", "maximum"], "DesignApproach"),
        (["$defs", "structElement", "properties", "kind", "maximum"], "StructKind"),
        (["$defs", "load", "properties", "kind", "maximum"], "LoadKind"),
        (["$defs", "soilPolygon", "properties", "edge_bc", "items", "maximum"], "BCType"),
        (["$defs", "soilPolygon", "properties", "edge_flow", "items", "maximum"], "FlowBCType"),
    ]
    for path, enum_name in expectations:
        if enum_name not in sizes:
            problems.append(f"project.hpp: cannot find enum class {enum_name}")
            continue
        want = sizes[enum_name] - 1
        got = bound(path)
        if got != want:
            problems.append(
                f"docs/k2d.schema.json: $defs.{'.'.join(path)} is {got}, but enum {enum_name} "
                f"has {sizes[enum_name]} values so the maximum must be {want} "
                f"(append-only enums widen the schema in the same commit)")

    if problems:
        for p in problems:
            print(f"  {p}")
        print(f"check_k2d_spec: {len(problems)} problem(s)")
        return 1
    print(f"check_k2d_spec: OK - {len(written)} keys locked across writer/reader/spec/schema; "
          f".k2d v{k2d_code}, .res v{res_code}; {len(expectations)} enum bounds verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
