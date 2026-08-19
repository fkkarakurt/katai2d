# The KATAI 2D file formats: `.k2d` and `.res`

This document is the written half of the input contract. The executable halves are the reader
(`kernel/io/include/katai/io/project_io.hpp`), the validator (`kernel/io/include/katai/io/validate.hpp`)
and the published JSON Schema (`docs/k2d.schema.json`). The CTest gate `check_k2d_spec` keeps all four
in agreement: every key the writer emits must appear in the field inventory below and in the JSON
Schema, every key documented here must still exist in the code, and the enum bounds in the schema must
equal the enums in `kernel/model/include/katai/model/project.hpp`. A hand-maintained format document
that can drift from the code would be a silent-wrong of its own kind; this one cannot drift silently.

Current `.k2d` version: **14** · Current `.res` version: **6**

Version history: **v2** adds line prescribed displacements (`disps` and the phase `disp`
activity flags). **v3** adds the anchor lock-off force (`anchors[i].prestress`), **v4** the
per-phase water conditions (`phases[i].water_override`, `wx`, `wy`), **v5**
the dilatancy cut-off (`materials[i].dilatancy_cutoff`, `e_max`), **v6** the prescribed
boundary flux (`polygons[i].edge_flux`, and `3` in `edge_flow`) and **v7** the per-phase
numerical controls (`tol`, `loadsteps`, `maxiter`), so that a published run carries the
numerics it was computed with, **v8** the staged-construction target and the undrained
switch (`mstage`, `ignoreund`) and **v9** the per-material undrained stiffness
(`materials[i].und_mode`, `nu_u`, `skempton_B`), which an older build would replace with
PLAXIS's default of 0.495 for every undrained material in the file, and **v10** the
Undrained (C) drainage type (`materials[i].drainage` = 4), which an older build would load
"for display as Drained" and solve with undrained parameters read as effective ones, and **v11**
wells and drains (`hydros`, `phases[i].hydro`), which an older build would ignore and solve the
same ground with no dewatering in it, and **v12** the cross permeability of walls and interfaces
(`structs[i].flow_barrier`, `hyd_res`), which an older build would ignore and compute the flow
straight through a cut-off wall, and **v13** the embedded beam's connection point
(`structs[i].conn`), which an older build would ignore and leave every pile top free of the soil
— the reading in which a point load at a pile head is carried by the ground beside it rather than
by the pile, and the same drawing settles several times as much, and **v14** the small-strain
history reset (`phases[i].resetsmall`), which an older build would ignore and run the phase on
the stiffness the earlier phases had already degraded, so the settlements and deflections it
reports are systematically larger than the ones the file asks for. Every bump
is deliberate and for the same reason: an older build reading the newer file would silently
drop the input and solve a *different* problem -- a wall with slack anchors deflects far more
than one that was tensioned against it -- so it must refuse the file instead.

## 1. The `.k2d` project file

A `.k2d` file is one JSON object, UTF-8 encoded, holding the complete project: geometry, material
data sets, structures, loads, water, boundary conditions and the staged-construction phase list.
Everything a run needs is in the file — an accelerogram, for example, is stored inline rather than
referenced, so a result is reproducible from the file alone.

### 1.1 Container rules

- **One strict JSON object.** The reader is a strict recursive-descent parser: it accepts objects,
  arrays, strings, numbers, `true`/`false`/`null`, rejects anything after the closing brace
  ("trailing characters") and rejects a root that is not an object.
- **Version guard.** The root key `katai2d` carries the format version. A file with `katai2d`
  missing, below 1 or above the version this build knows is refused: *"not a KATAI 2D project file
  (or a newer version)"*.
- **Missing key = documented default.** Every key except `katai2d` is optional. An absent key loads
  the default listed in the inventory below. This is the backward-compatibility rule: files written
  before a field existed keep loading, and the default reproduces the old behaviour.
- **Unknown key = ignored.** The reader looks keys up by name and never rejects a key it does not
  know. This is the forward-compatibility rule: a file from a newer minor revision still loads,
  minus what this build cannot represent.
- **Enum values are file-stable.** Enumerations are stored as integers and new values are only ever
  APPENDED. An integer outside the known range is clamped to a safe display value (soil model →
  `0` Linear elastic, drainage → `0` Drained) and **reported** through the reader's notes channel —
  a front end that intends to solve must treat such a note as an error. Enums that are not clamped
  at read time (phase type, structure kind, load kind, waveform, design approach) are caught by the
  validator instead.
- **Doubles round-trip exactly.** The writer emits numbers with `%.17g`, which reproduces IEEE
  doubles bit-exactly. NaN and infinities never enter a file: the writer maps them to `0`.
- **Duplicate keys: the first occurrence wins.** The writer never emits duplicates.
- **String escapes.** The writer escapes `"`, `\`, and control characters (`\n`, `\t`, `\r`,
  `\uXXXX` for the rest) and emits all other bytes raw (UTF-8). The reader additionally accepts
  `\uXXXX` escapes and decodes the Latin-1 range; a code point above `ÿ` is replaced with `?`.
- **Version bump policy.** Adding a field or appending an enum value does NOT bump `katai2d` — the
  missing-key and append-only rules cover both directions. The version bumps only when an existing
  key changes meaning, unit or shape, which is a breaking change and should be rare and deliberate.

### 1.2 Validation

Structure (types, array shapes) is described by `docs/k2d.schema.json`. Everything beyond structure
— index cross-references, parameter bounds, silently-ignored input — is the job of
`katai::io::validate_project()`, which reports every violation at its exact field path
(`materials[2].Rf`) in engineer-readable language with the offending value echoed. The JSON Schema
deliberately does **not** duplicate the validator's parameter bounds: one source of truth per rule.
Model × drainage × phase feasibility (e.g. Soft Soil Creep under a Consolidation phase) is likewise
not duplicated anywhere: the engine refuses those combinations at solve time with its own tested
messages.

### 1.3 Field inventory

Types: `int`, `num` (double), `bool`, `str`, `num[]` (array of numbers), `int[]`, `01[]` (array of
0/1 activity flags). "Default" is what an ABSENT key loads as — the reader's value, which in three
flagged places differs from the in-memory default of a freshly created object.

#### Root object

<!-- field-inventory:begin -->

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `katai2d` | int | — | *required* | Format version; this build writes 12 |
| `name` | str | — | `"Untitled project"` | Project display name |
| `axisymmetric` | bool | — | `false` | `false` = plane strain; `true` = axisymmetric (x is the radius, x = 0 the symmetry axis) |
| `initial_procedure` | int | — | `0` | How the initial phase establishes the in-situ state: 0 K0 procedure, 1 Gravity loading, 2 Safety (φ-c reduction of the gravity state; intended for single-phase runs) |
| `x_min` | num | m | `0` | Model extent, left edge |
| `x_max` | num | m | `40` | Model extent, right edge |
| `y_min` | num | m | `0` | Model extent, base |
| `y_max` | num | m | `20` | Model extent, top |
| `has_water` | bool | — | `false` | Groundwater present. NOTE: an absent key loads `false`, although a freshly created project defaults to `true` — the reader's value is the contract |
| `wx` | num[] | m | `[]` | Water polyline x-coordinates (with `wy`; needs ≥ 2 points to take effect) |
| `wy` | num[] | m | `[]` | Water polyline y-coordinates (elevation head) |
| `mesh` | object | — | defaults | Global mesh-generation settings (object below); the file determines its own discretization |
| `materials` | array | — | `[]` | Soil material data sets (objects below) |
| `plates` | array | — | `[]` | Plate material data sets |
| `anchors` | array | — | `[]` | Anchor material data sets |
| `geogrids` | array | — | `[]` | Geogrid material data sets |
| `embedded` | array | — | `[]` | Embedded-beam (pile row) material data sets |
| `polygons` | array | — | `[]` | Soil regions (user-drawn polygons) |
| `structs` | array | — | `[]` | Structural elements (lines) |
| `loads` | array | — | `[]` | External loads |
| `disps` | array | — | `[]` | Line prescribed displacements (objects below) |
| `hydros` | array | — | `[]` | Wells and drains — hydraulic conditions inside the model (objects below) |
| `initial` | object | — | defaults | The initial phase (activation flags; its `type` is ignored — `initial_procedure` selects how the initial state is established) |
| `phases` | array | — | `[]` | Staged-construction phases after the initial phase, in run order |

#### Mesh object (`mesh`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `elem_size` | num | m | `2` | Target element edge length; the mesher's target triangle area is 0.5 · elem_size² |
| `order` | int | — | `6` | Element type: 6 = quadratic tri6, 15 = quartic tri15 |
| `auto_refine` | bool | — | `true` | Automatic ~2× refinement near structural elements and loads |

#### Soil material object (`materials[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `name` | str | — | `"New material"` | Data set name |
| `model` | int | — | `0` | Soil model: 0 Linear elastic, 1 Mohr-Coulomb, 2 Hardening Soil, 3 HS small, 4 Soft Soil, 5 Soft Soil Creep. NOTE: an absent key loads 0, although a freshly created material defaults to 1 |
| `drainage` | int | — | `0` | 0 Drained, 1 Undrained (A), 2 Non-porous, 3 Undrained (B), 4 Undrained (C) — the last is a **total stress** analysis: `E`/`nu` are the undrained pair, `c` is s_u with φ = 0, and no pore pressure is generated or carried (PLAXIS MMM §2.7) |
| `color` | num[] | — | `[0.85,0.78,0.55]` | Display colour, RGB in 0..1 (3 values) |
| `gamma_unsat` | num | kN/m³ | `17` | Unit weight above the phreatic level |
| `gamma_sat` | num | kN/m³ | `20` | Unit weight below the phreatic level |
| `e_init` | num | — | `0.5` | Initial void ratio |
| `E` | num | kN/m² | `1.3e4` | Young's modulus (Linear elastic / Mohr-Coulomb) |
| `nu` | num | — | `0.3` | Poisson's ratio (Linear elastic / Mohr-Coulomb) |
| `c` | num | kN/m² | `1` | Effective cohesion c′ (Undrained (B): the undrained strength su) |
| `phi` | num | ° | `30` | Effective friction angle φ′ |
| `psi` | num | ° | `0` | Dilatancy angle ψ |
| `E_inc` | num | kN/m²/m | `0` | Stiffness increase per metre depth below `y_ref` |
| `c_inc` | num | kN/m²/m | `0` | Cohesion increase per metre depth below `y_ref` |
| `y_ref` | num | m | `0` | Reference level for the increments |
| `tension_cutoff` | bool | — | `true` | Rankine tension cut-off active |
| `dilatancy_cutoff` | bool | — | `false` | Stop dilatancy at the critical void ratio (PLAXIS MMM Eq. 5.16b) |
| `e_max` | num | — | `1` | Critical (maximum) void ratio; read when `dilatancy_cutoff` is set |
| `tensile_strength` | num | kN/m² | `0` | Allowed tensile strength σ_t |
| `E50ref` | num | kN/m² | `3e4` | HS: secant stiffness at p_ref |
| `Eoedref` | num | kN/m² | `3e4` | HS: oedometer stiffness at p_ref |
| `Eurref` | num | kN/m² | `9e4` | HS: unload-reload stiffness at p_ref |
| `m` | num | — | `0.5` | HS: stress-dependency power |
| `nu_ur` | num | — | `0.2` | HS / Soft Soil: unload-reload Poisson's ratio |
| `p_ref` | num | kN/m² | `100` | HS: reference pressure |
| `Rf` | num | — | `0.9` | HS: failure ratio q_f/q_a |
| `k0nc_auto` | bool | — | `true` | K0ⁿᶜ from 1 − sin φ′ |
| `k0nc` | num | — | `0.5` | K0ⁿᶜ when `k0nc_auto` is false |
| `G0ref` | num | kN/m² | `1.2e5` | HS small: small-strain shear modulus at p_ref |
| `gamma07` | num | — | `1.5e-4` | HS small: shear strain at 0.722 G0 |
| `lamstar` | num | — | `0.10` | Soft Soil: modified compression index λ* |
| `kapstar` | num | — | `0.02` | Soft Soil: modified swelling index κ* |
| `mustar` | num | — | `0.005` | Soft Soil Creep: modified creep index μ* |
| `kx` | num | m/day | `1` | Horizontal permeability |
| `ky` | num | m/day | `1` | Vertical permeability |
| `und_mode` | int | — | `0` | How the pore fluid's stiffness is defined for Undrained (A)/(B): 0 = `nu_u` entered, 1 = Skempton's `skempton_B` (PLAXIS MMM §2.4) |
| `nu_u` | num | — | `0.495` | Equivalent undrained Poisson ratio (read when `und_mode` = 0); must satisfy ν' < ν_u < 0.5 |
| `skempton_B` | num | — | `0` | Skempton's B (read when `und_mode` = 1); must lie in (0, 1) |
| `gw_ga` | num | 1/m | `14.5` | van Genuchten g_a (inverse air-entry) |
| `gw_gn` | num | — | `2.68` | van Genuchten g_n (> 1) |
| `gw_gl` | num | — | `0.5` | Mualem pore-connectivity g_l |
| `gw_Sres` | num | — | `0.10` | Residual saturation |
| `rinter_rigid` | bool | — | `true` | Interface strength rigid (R_inter = 1) |
| `Rinter` | num | — | `1` | Interface strength factor when not rigid |
| `k0_auto` | bool | — | `true` | K0 from 1 − sin φ′ |
| `k0` | num | — | `0.5` | K0 when `k0_auto` is false |
| `oc_mode` | int | — | `0` | Stress history: 0 none, 1 OCR, 2 POP |
| `OCR` | num | — | `1` | Overconsolidation ratio (`oc_mode` = 1) |
| `POP` | num | kN/m² | `0` | Pre-overburden pressure (`oc_mode` = 2) |

#### Plate material object (`plates[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `elastoplastic` | bool | — | `false` | Check the plastic capacities Mp/Np |
| `EA` | num | kN/m | `5e6` | Axial stiffness |
| `EI` | num | kN·m²/m | `8.5e3` | Bending stiffness |
| `w` | num | kN/m/m | `0` | Plate weight |
| `Mp` | num | kN·m/m | `0` | Plastic moment capacity (0 = elastic) |
| `Np` | num | kN/m | `0` | Plastic axial capacity (0 = elastic) |

with `name` (default `"Plate"`), `color` and `nu` (default `0`) as in the soil material table.

#### Anchor material object (`anchors[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `Fmax_tens` | num | kN | `0` | Tension capacity (0 = unlimited) |
| `Fmax_comp` | num | kN | `0` | Compression capacity (0 = unlimited) |
| `Lspacing` | num | m | `1` | Out-of-plane spacing |
| `prestress` | num | kN | `0` | Lock-off force of one anchor, tension-positive (0 = installed slack) |

with `name` (default `"Anchor"`), `color`, `elastoplastic` and `EA` (default `1e5` kN) as above.

`prestress` is applied when the anchor is first active and the anchor is an elastic spring from
that state afterwards, so its force follows the wall rather than staying at the lock-off value.
The stated force is per anchor; the analysis divides it by `Lspacing`, exactly as it does `EA`
and the capacities.

#### Geogrid material object (`geogrids[i]`)

`name` (default `"Geogrid"`), `color`, `elastoplastic`, `EA` (default `1e3` kN/m) and `Np`
(tension cap, kN/m, default `0` = unlimited) as above.

#### Embedded-beam material object (`embedded[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `gamma` | num | kN/m³ | `24` | Beam unit weight |
| `diameter` | num | m | `0.40` | Pile diameter (circular) |
| `Tskin_max` | num | kN/m | `0` | Skin resistance cap (0 = unlimited) |
| `Fmax_base` | num | kN | `0` | Base resistance (0 = unlimited) |

with `name` (default `"Embedded beam"`), `color`, `E` (default `3e7` kN/m²) and `Lspacing`
(default `2.5` m) as above.

#### Soil region object (`polygons[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `material` | int | — | `-1` | Index into `materials` (−1 = unassigned; the validator rejects it) |
| `coarseness` | num | — | `1` | Local mesh density factor (1 = global, 0.5 = twice as fine) |
| `x` | num[] | m | `[]` | Polygon vertex x-coordinates (implicitly closed; ≥ 3) |
| `y` | num[] | m | `[]` | Polygon vertex y-coordinates |
| `edge_bc` | int[] | — | `[]` | Per-edge deformation BC (edge j runs vertex j → j+1): 0 Free, 1 Normally fixed, 2 Horizontally fixed, 3 Vertically fixed, 4 Fully fixed. Empty or one per vertex |
| `edge_flow` | int[] | — | `[]` | Per-edge flow BC: 0 Closed, 1 Prescribed head, 2 Seepage face, 3 Prescribed flux |
| `edge_head` | num[] | m | `[]` | Prescribed head per edge (read where `edge_flow` is 1) |
| `edge_flux` | num[] | m/day | `[]` | Prescribed boundary flux per edge, inflow positive (read where `edge_flow` is 3) |

with `name` (default `"Soil"`) as above.

#### Structural element object (`structs[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `kind` | int | — | `0` | 0 Plate, 1 Anchor, 2 Geogrid, 3 Embedded beam, 4 Interface |
| `x1` | num | m | `0` | Start point x |
| `y1` | num | m | `0` | Start point y |
| `x2` | num | m | `0` | End point x |
| `y2` | num | m | `0` | End point y |
| `iface_pos` | bool | — | `false` | Positive-side interface attached |
| `iface_neg` | bool | — | `false` | Negative-side interface attached |
| `iface_material` | int | — | `-1` | Soil material override for the interfaces (−1 = adjacent soil) |
| `conn` | int | — | `0` | **Embedded beam only** — how the connection point (the pile top: the end with the highest y, or for an exactly horizontal pile the lowest x) is attached (PLAXIS Reference §5.6.3). `0` hinged: the beam's top translations *are* the soil's there, "the same displacement, but not necessarily the same rotation" — PLAXIS's default when no structure shares the point, and what makes a pile loadable at its head. `1` free: coupled to the soil through the skin springs only, which is what PLAXIS sets for a grout body so a ground anchor does not shed axial force at the connection. The mesher carries the connection point as a node either way, so switching this changes the physics and not the mesh |
| `flow_barrier` | int | — | `0` | Cross permeability in a groundwater calculation: 0 fully permeable (no effect on flow), 1 impermeable (the two sides get separate pore-pressure DOFs), 2 semi-permeable (PLAXIS Ref Table 5-2). Plates and interfaces only |
| `hyd_res` | num | day | `0` | Hydraulic resistance d/k of a semi-permeable barrier: the head difference divided by the discharge per unit area of wall |

with `name` (default `"Element"`), `material` (index into the kind's material list, default `-1`)
and `coarseness` (default `1`) as above. An Interface element takes its strength from the adjacent
soil (or `iface_material`); its `material` field is ignored.

#### Load object (`loads[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `qx1` | num | kN or kN/m | `0` | x-component at the point / start of the line |
| `qy1` | num | kN or kN/m | `0` | y-component at the point / start. NOTE: an absent key loads 0, although a freshly created load defaults to −10 |
| `qx2` | num | kN/m | `0` | x-component at the end (Distributed; varies linearly) |
| `qy2` | num | kN/m | `0` | y-component at the end. Same reader-default note as `qy1` |

with `kind` (0 Point at (x1, y1) in kN; 1 Distributed along (x1, y1)–(x2, y2) in kN/m, default
`0`), `name` (default `"Load"`), `x1`/`y1`/`x2`/`y2` and `coarseness` as above.

#### Prescribed-displacement object (`disps[i]`)

Every mesh node on the line (x1, y1)–(x2, y2) gets the *set* components imposed, ramped
0 → value over the phase like a load; unset components stay free. A set component with
value 0 is a rigid support line. Activated per phase through the phase `disp` flags;
static (Plastic) phases only in this build. The classic use is a displacement-controlled
rigid footing whose force is read from the reactions.

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `set_ux` | bool | — | `false` | The horizontal component is prescribed |
| `ux` | num | m | `0` | Horizontal displacement (when `set_ux`) |
| `set_uy` | bool | — | `true` | The vertical component is prescribed |
| `uy` | num | m | `0` | Vertical displacement (when `set_uy`) |

with `name` (default `"Displacement"`), `x1`/`y1`/`x2`/`y2` and `coarseness` as above.

#### Hydraulic-condition object (`hydros[i]`)

A well or a drain drawn *inside* the model (PLAXIS Reference §5.9). Its line is embedded in the
mesh like a load line. A **well** prescribes a discharge spread along it and stops extracting once
the head reaches `h_min`; a **drain** holds the head at `head` — a *normal* drain only takes water
away, a *vacuum* drain holds its head in both directions. Activated per phase through the phase
`hydro` flags. Read by the groundwater-flow calculation; in consolidation and fully-coupled phases
a drain sets the excess pore pressure to zero and a well is reported as not applied (`K2D-A009`).

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `kind` | int | — | `0` | 0 Well (a discharge), 1 Drain (a head) |
| `behaviour` | int | — | `0` | Well: 0 Extraction, 1 Infiltration. Drain: 0 Normal, 1 Vacuum |
| `q` | num | m³/day/m | `0` | Well discharge \|Q\| per metre out of plane (> 0) |
| `h_min` | num | m | `0` | Well: the lowest head it can draw the ground down to |
| `head` | num | m | `0` | Drain: the head it holds |

with `name` (default `"Well"`), `x1`/`y1`/`x2`/`y2` and `coarseness` as above.

#### Phase object (`initial` and `phases[i]`)

| Key | Type | Unit | Default | Meaning |
|---|---|---|---|---|
| `type` | int | — | `0` | 0 Plastic, 1 Safety, 2 Consolidation, 3 Transient flow, 4 Fully coupled, 5 Dynamic |
| `duration` | num | day (s for Dynamic) | `1` | Time interval of a time-dependent phase; creep time of a chained Plastic phase |
| `steps` | int | — | `25` | Time steps over `duration` (Newmark steps for Dynamic) |
| `design` | int | — | `0` | Design approach: 0 None, 1 EC7 DA1-C1, 2 EC7 DA1-C2, 3 EC7 DA2, 4 EC7 DA3, 5 TBDY 2018 static, 6 TBDY 2018 seismic |
| `seiswave` | int | — | `0` | Dynamic input waveform: 0 Harmonic, 1 Ricker, 2 Record |
| `seisamp` | num | m/s² | `1` | Peak base acceleration (Record: scale factor, 1 = as recorded) |
| `seisfreq` | num | Hz | `2` | Dominant frequency (Harmonic/Ricker) |
| `damp` | num | — | `0.05` | Rayleigh damping ratio at the two target frequencies |
| `rayf1` | num | Hz | `1` | First Rayleigh target frequency |
| `rayf2` | num | Hz | `8` | Second Rayleigh target frequency (> `rayf1`) |
| `seisff` | num 0/1 | — | `0` | Free-field lateral boundaries (Lysmer) — stored as 0/1, not JSON booleans |
| `dynnl` | num 0/1 | — | `0` | Fully nonlinear Newmark solve (plastification during shaking) |
| `seiscb` | num 0/1 | — | `0` | Compliant (absorbing) base, input as the upward wave |
| `ec8on` | num 0/1 | — | `0` | EC8 design-spectrum overlay active |
| `ec8agr` | num | g | `0.2` | EC8 reference peak ground acceleration a_gR |
| `ec8gi` | num | — | `1` | EC8 importance factor γ_I |
| `ec8gnd` | int | — | `2` | EC8 ground type: 0 A … 4 E |
| `ec8typ` | int | — | `0` | EC8 spectrum: 0 Type 1, 1 Type 2 |
| `rec` | num[] | m/s² | `[]` | Inline accelerogram, uniformly sampled (written only when non-empty) |
| `recdt` | num | s | `0.02` | Accelerogram sampling interval |
| `tbdyss` | num | — | `1` | TBDY 2018 short-period spectral coefficient S_S |
| `tbdys1` | num | — | `0.4` | TBDY 2018 1-second spectral coefficient S_1 |
| `siteclass` | int | — | `2` | TBDY site class: 0 ZA … 4 ZE |
| `poly` | 01[] | — | `[]` | Soil-region activity flags by index; a MISSING entry counts as active |
| `struct` | 01[] | — | `[]` | Structural-element activity flags |
| `load` | 01[] | — | `[]` | Load activity flags |
| `disp` | 01[] | — | `[]` | Prescribed-displacement activity flags |
| `hydro` | 01[] | — | `[]` | Well / drain activity flags |
| `water_override` | bool | — | `false` | Use this phase's own phreatic line instead of the project's |
| `wx` | num[] | m | `[]` | Phase phreatic polyline, x (used when `water_override`) |
| `wy` | num[] | m | `[]` | Phase phreatic polyline, y (used when `water_override`) |
| `mstage` | num | — | `1` | Fraction of this phase's staged change to apply (PLAXIS Σ Mstage). Staged (non-initial) phases only. Written only when below 1 |
| `ignoreund` | bool | — | `false` | Solve Undrained (A)/(B) materials as drained in this phase (PLAXIS "Ignore und. behaviour"); strength parameters unchanged. Written only when true |
| `resetsmall` | bool | — | `false` | Clear the Hardening Soil small-strain history at the start of this phase (PLAXIS "Reset small strain", Material Models Manual sec. 7.6), so the soil meets the phase at G0 instead of the stiffness earlier phases degraded it to. Stress, shear hardening `gamma_p` and the preconsolidation pressure are carried over untouched — the option resets the *history*, not the state. Staged (non-initial) phases only; raises `K2D-M005` saying how many stress points it cleared, or that no small-strain material was present to clear. Written only when true |
| `tol` | num | — | *by material class* | Tolerated relative force residual (PLAXIS "Tolerated error"). Written only when set |
| `loadsteps` | int | — | *by material class* | Load increments for this phase. Not PLAXIS's "Max steps": KATAI splits the load into a fixed number of increments (with adaptive cut-back), it does not step automatically to a cap. Written only when set |
| `maxiter` | int | — | *by phase strategy* | Newton iterations per increment (PLAXIS "Max iterations"). Written only when set |

with `name` (default `"Phase"`) as above.

<!-- field-inventory:end -->

## 2. The `.res` results file

Binary on purpose: FEM results are large numeric arrays, and JSON would cost ~20 bytes per double
plus a re-parse. Layout (little-endian x86, doubles stored raw):

```
"K2DR" magic · u32 version · u64 model_hash · mesh (stored once — phases share it) ·
u32 phase count · per phase: flags + scalars + message + displacements/stresses/pore/activity
+ structural force diagrams + interface results
```

- `model_hash` is FNV-1a over the CANONICAL project JSON (`project_to_json` of the loaded model),
  so results are refused as **stale** when the model no longer matches them — even after
  formatting-only edits to the project file.
- Version history (a reader accepts 1..current; absent fields read as `false`/`0`, which is correct
  for files that predate the feature): **v2** interface results (τ/σ_n/slip) · **v3** seismic
  envelope flags · **v4** superposed design action + Coulomb utilisation · **v5**
  `slip_checked` (the envelope came from the nonlinear Coulomb branch) · **v6** `stopped_by`,
  the reason a solve stopped short of full load, without which a reopened result cannot say whether
  its load factor is a capacity or the point where the iteration budget ran out.
- Committed Gauss states are NOT stored: restored results are for viewing and post-processing;
  continuing a staged run re-calculates. This is honest and keeps the file an order of magnitude
  smaller.
- A corrupt or truncated file is refused with a message, never partially loaded; array sizes are
  checked against the stored mesh.

Round-trip and refusal behaviour are pinned by `tests/test_results_io.cpp`; the project round trip
by `tests/test_project_io.cpp`; the validator rules by `tests/test_k2d_validator.cpp`.
