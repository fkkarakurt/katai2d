"""Stub completeness: every public `_core` name appears in the generated stub.

_core.pyi is derived from the built extension by nanobind's stubgen (a
POST_BUILD step), so it cannot drift from the bindings by construction.
What CAN silently break is the generation itself -- a stubgen upgrade that
starts skipping a construct, or the step being dropped from the build.
This test pins the contract the stub exists for: an IDE user must be able
to complete every public name the module actually exports.
"""

import re
import sys
from pathlib import Path

import katai._core as core

stub = Path(core.__file__).with_name("_core.pyi")
assert stub.exists(), f"missing {stub} -- the stubgen POST_BUILD step did not run"
text = stub.read_text(encoding="utf-8")

# DEFINED top-level names, not substrings: a name that merely occurs in some
# docstring must not count as covered, so the check parses definition lines.
defined = set(re.findall(r"^(?:class|def)\s+(\w+)", text, re.M))
defined |= set(re.findall(r"^(\w+)\s*[:=]", text, re.M))

public = [n for n in dir(core) if not n.startswith("_")]
missing = [n for n in public if n not in defined]
if missing:
    for name in missing:
        print(f"  _core.{name} is exported but not defined in _core.pyi")
    print(f"stub_check: {len(missing)} name(s) missing")
    sys.exit(1)

print(f"ok: _core.pyi defines all {len(public)} public names")
