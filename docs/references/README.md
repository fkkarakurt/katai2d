# Reference library

The primary sources this implementation rests on. The purpose is that an implementation or
verification question is answered *here* rather than by searching again: a formulation whose source
cannot be named is not a formulation, and a secondary source — a blog, a forum, a summary — is never
the basis of one.

## Backbone sources

| Short | Source | Used for |
|---|---|---|
| **PLAXIS-MM** | PLAXIS 2D 2025.1 *Material Models Manual* (Seequent, free PDF) | matching PLAXIS behaviour exactly: Mohr-Coulomb, Hardening Soil, HS-small, Soft Soil, Soft Soil Creep, Undrained A/B/C, interfaces. The primary external reference. |
| **PLAXIS-Sci** | PLAXIS 2D 2025.1 *Scientific Manual* | element formulations, integration, dynamics and absorbing boundaries |
| **PLAXIS-Val** | PLAXIS 2D 2025.1 *Validation Manual* | published benchmark cases with values to reproduce |
| **P&Z-1 / P&Z-2** | Potts & Zdravković, *Finite Element Analysis in Geotechnical Engineering* — Vol 1 Theory (1999), Vol 2 Application (2001), Thomas Telford | the standard reference for geotechnical finite elements: theory and application |
| **SGM** | Smith, Griffiths & Margetts, *Programming the Finite Element Method*, 5th ed. (2014), Wiley | working formulations for elasticity, plasticity, consolidation and flow |
| **Simo-Hughes** | Simo & Hughes, *Computational Inelasticity* (1998), Springer | return mapping and the integration of plasticity |
| **dSN** | de Souza Neto, Perić & Owen, *Computational Methods for Plasticity* (2008), Wiley | Mohr-Coulomb return mapping and practical algorithms |
| **ZT** | Zienkiewicz & Taylor, *The Finite Element Method* (6th–7th ed.) | general finite element theory, element formulation, numerical integration |

Official location of the PLAXIS manuals (version in the path may advance):
`https://files.seequent.com/PLAXIS/Manuals/PLAXIS_2D/English/`

## Formulation notes

One note per model or analysis type. Each states the governing equations, the source and clause they
come from, and the verification that pins the implementation to them.

| Area | Notes |
|---|---|
| Constitutive | [constitutive-models.md](constitutive-models.md) · [mohr-coulomb-formulation.md](mohr-coulomb-formulation.md) · [hardening-soil-formulation.md](hardening-soil-formulation.md) · [hssmall-formulation.md](hssmall-formulation.md) · [soft-soil-formulation.md](soft-soil-formulation.md) · [soft-soil-creep-formulation.md](soft-soil-creep-formulation.md) · [material-model-architecture.md](material-model-architecture.md) |
| Stress state | [effective-stress-formulation.md](effective-stress-formulation.md) · [initial-stress-k0.md](initial-stress-k0.md) |
| Elements and meshing | [elements-and-meshing.md](elements-and-meshing.md) · [meshing-design.md](meshing-design.md) · [mesh-sizing.md](mesh-sizing.md) |
| Structural | [structural-plate-formulation.md](structural-plate-formulation.md) · [interface-formulation.md](interface-formulation.md) · [embedded-beam-formulation.md](embedded-beam-formulation.md) |
| Water | [seepage-formulation.md](seepage-formulation.md) · [transient-unsaturated-flow-formulation.md](transient-unsaturated-flow-formulation.md) · [consolidation-formulation.md](consolidation-formulation.md) |
| Analysis | [analysis-and-structural.md](analysis-and-structural.md) · [dynamic-seismic-formulation.md](dynamic-seismic-formulation.md) |
| Design codes | [design-codes-ec7-tbdy.md](design-codes-ec7-tbdy.md) · [tbdy-2018-seismic.md](tbdy-2018-seismic.md) |
| Application | [gui-design.md](gui-design.md) |
| Verification | [validation-benchmarks.md](validation-benchmarks.md) · [multi-code-validation-plan.md](multi-code-validation-plan.md) · [literature-review.md](literature-review.md) · [plaxis-gap-analysis.md](plaxis-gap-analysis.md) |

## Use of sources

These are references for equations and algorithms. No third-party code is copied: the SGM listings are
Fortran and licensed with the book, so they are read and the C++ here is written independently. The
PLAXIS manuals document formulations, and the implementation rests on the published equations, not on
any part of that product.
