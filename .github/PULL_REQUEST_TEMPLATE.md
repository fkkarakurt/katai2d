# What this changes, and why

<!-- One logical change per pull request. -->

# Verification

<!-- Required when the change touches physics. Which oracle backs it: a
     closed-form solution (stated in full in the test), an independent
     computation path that shares no code with the solver, or a published
     benchmark (primary source cited)? Self-consistency does not count. -->

# Checklist

- [ ] The full suite passes (`.\scripts\build.ps1 -Configure -Test`) — the suite, not a subset, is the specification.
- [ ] A physics change carries a test with its reference declaration, and the matrix regenerates cleanly (`python scripts/check_references.py --write`).
- [ ] A `.k2d` format change appears in all four places: writer, reader, `docs/k2d-format.md`, `docs/k2d.schema.json`.
- [ ] English only in code, comments, messages and documentation.
- [ ] Commits are signed off (`git commit -s`, Developer Certificate of Origin).
