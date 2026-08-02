# KATAI 2D — documentation

A 2D finite element engine for geotechnical engineering: soil, structural elements and groundwater
in one staged-construction workflow, driven from the command line or from Python.

> **Working principle: simplest to most complex.** Every step is verified against an independent
> oracle before the next one starts, and a capability that is not verified is not in the product.
> A plausible wrong number is treated as worse than a crash, because it does not announce itself.

## Map

| Document | Contents |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | the layer contract, module tree, include and namespace conventions, and an explicit list of the places where the tree does not yet match the design |
| [DECISIONS.md](DECISIONS.md) | locked technical decisions and open questions |
| [ENVIRONMENT.md](ENVIRONMENT.md) | toolchain, build script, dependency handling |
| [how-to-design-codes.md](how-to-design-codes.md) | how a design-code catalogue is added |
| [references/](references/) | formulation notes: one document per model or analysis type, each tied to its primary source |
| [validation/](validation/) | the verification and validation record, including the cross-code comparison table and the performance baseline |

## Verification

The suite is the specification. Every claim in `validation/` is backed by a test that runs in CTest,
and every test is checked against a closed-form solution, an independent path that shares no code
with the solver, or an external code or publication. Self-consistency does not count as
verification.

```
python scripts/check_language.py --report      # non-English prose census
python scripts/check_architecture.py --list    # module graph resolved from the include statements
ctest --test-dir build/<preset>                # the full suite
```
