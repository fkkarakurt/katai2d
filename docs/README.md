# KATAI 2D — documentation

A 2D finite element engine for geotechnical engineering: soil, structural elements and groundwater
in one staged-construction workflow, driven from the command line or from Python.

> **Working principle: simplest to most complex.** Every step is verified against an independent
> oracle before the next one starts, and a capability that is not verified is not in the product.
> A plausible wrong number is treated as worse than a crash, because it does not announce itself.

## Map

| Document | Contents |
|---|---|
| [k2d-format.md](k2d-format.md) | the `.k2d` input format: field inventory, units, versioning — the one contract the file, the script and the solver share |
| [k2d.schema.json](k2d.schema.json) | the same contract, machine-readable (JSON Schema) |
| [validation/verification-matrix.md](validation/verification-matrix.md) | the verification matrix, generated from the reference declarations inside the tests |
| [validation/references.bib](validation/references.bib) | every primary source the matrix cites, generated alongside it |
| [diagnostics.md](diagnostics.md) | every `K2D-*` code the engine can raise, what it means and what to do about it |
| [validation/three-published-benchmarks.md](validation/three-published-benchmarks.md) | three published benchmarks reproduced letter-for-letter and solved with the `katai` CLI, published value beside computed value |
| [validation/plaxis-2d-validation-comparison.md](validation/plaxis-2d-validation-comparison.md) | four cases from a published validation manual, rebuilt from their problem statements |
| [validation/numerical-uncertainty.md](validation/numerical-uncertainty.md) | how far a computed number sits from the exact solution of its own equations, measured rather than asserted |

The spec, the schema and the matrix are pinned to the code by the test suite (`check_k2d_spec`,
`test_reference_registry`): they cannot drift from what the program actually reads, writes and
verifies. The longer-form formulation notes and validation records are working documents and live
outside the repository; they are being prepared for the project website.

## Verification

The suite is the specification. Every claim in the verification matrix is backed by a test that
runs in CTest, and every test is checked against a closed-form solution, an independent path that
shares no code with the solver, or an external code or publication. Self-consistency does not
count as verification.

```
python scripts/check_language.py --report      # non-English prose census
python scripts/check_architecture.py --list    # module graph resolved from the include statements
ctest --test-dir build/<preset>                # the full suite
```
