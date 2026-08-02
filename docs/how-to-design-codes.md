# Running a design-code check (Eurocode 7 / TBDY 2018)

A code check in KATAI 2D is a characteristic analysis with partial factors applied. What follows is
how to run one.

## Two kinds of approach

Codes reduce capacity in two different ways, and KATAI treats them differently because they are not
the same operation:

- **Material-factored** (EC7 **DA1-C2**, **DA3**): the soil's strength is reduced (`c'/γ_c'`,
  `tanφ'/γ_φ'`) and the problem is solved again. This maps directly onto the finite element method,
  and the result is an **Over-Design Factor (ODF)**.
- **Resistance-factored** (EC7 **DA2**, **TBDY 2018**): the factor is applied to the computed
  **resistance** (`R_d = R_k/γ_Rv`, checked as `E_d ≤ R_d`). The finite element solution stays
  characteristic and the factors live in the reporting layer.

> Most common in practice: **DA3** for general stability, slopes and excavations; **DA2** or **TBDY**
> for bearing capacity.

## Step by step: stability check with EC7 DA3 (or DA1-C2)

1. Draw the geometry, materials and loads as for any model.
2. Switch to **staged construction** and add a phase with **Add phase**.
3. Set the phase **Type** to **Safety (phi-c reduction)**.
4. In the same phase, set **Design approach** to **Eurocode 7 - DA3** (or DA1 Comb. 2).
5. **Calculate**.

The reported number is then an **Over-Design Factor**:

- **ODF ≥ 1.0** — the design satisfies the ultimate limit state.
- **ODF < 1.0** — it does not: strengthen, change the geometry, or add support.

Why a Safety phase? Under DA3 the strength is reduced in advance by the M2 set (γ = 1.25), and the
Safety analysis then finds the *additional* factor available on top of that before collapse. That
additional factor is the ODF, so `ODF = FoS_characteristic / 1.25`.

## EC7 DA2 and TBDY 2018 (resistance-factored)

Here the finite element solution runs with characteristic values and the report carries the factors
for the `E_d ≤ R_d = R_k/γ_Rv` check. TBDY 2018 Table 16.2 gives bearing capacity **1.40**, sliding
**1.10** and passive resistance **1.40** — identical to EC7 DA2. Select the approach in the phase and
calculate; the factors appear in Output → Report.

## Report

After any phase with a design approach, **Output → Detailed results & report** lists, with the clause
cited: the approach selected, the reference (EN 1997-1 / TBDY 2018 §16.8.2), the partial factors
applied, and the pass or fail verdict on the ODF.

## Reference

Factor values and their sources: `references/design-codes-ec7-tbdy.md` (EN 1997-1 Annex A and
TBDY 2018 Chapter 16, official gazette of 18 March 2018).
