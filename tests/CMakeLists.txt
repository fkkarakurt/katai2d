# tests -- unit and regression. No framework: CTest plus assertions is enough, and
# every test is a plain executable that returns non-zero on failure.
#
# katai_add_test is defined in cmake/KataiAddTest.cmake, which the root includes whether or not
# this directory is added -- see the header there for what the function does and why it moved.
#
# KATAI_CORE_TESTS_MATCHING narrows THIS directory only. It is copied into a directory-scoped
# name so that a superbuild adding its own tests after this one is unaffected by it.
set(KATAI_TEST_FILTER "${KATAI_CORE_TESTS_MATCHING}")

katai_add_test(test_math_smoke)                    # Eigen smoke (P0.1)
katai_add_test(test_sparse_matrix)                 # CSR + builder (P0.2)
katai_add_test(test_grid_convergence)              # KV-NUM-004: the discretisation-error estimator itself (Roache 1994 / Celik et al. 2008 GCI) against manufactured power-law triplets whose exact answer is known
# Sparse direct solvers: every backend against a dense reference, and backend
# against backend. Always built, in every configuration, because it is the test
# that says a result does not depend on which solver the program contains.
katai_add_test(test_linsolve katai::linsolve)
katai_add_test(test_tri6 katai::fem)              # tri6 element + LE material (P0.3)
katai_add_test(test_tri15 katai::fem)             # tri15 quartic element kernels (P1.5)
katai_add_test(test_axisymmetric katai::fem)      # axisymmetric element kernels (P1.6)
katai_add_test(test_axisym_mc katai::materials)         # axisymmetric MC material point + tangent (P1.6)
katai_add_test(test_mesh katai::mesh katai::fem)  # structured tri6 mesh (P0.4+P0.5)
katai_add_test(test_predicates katai::mesh)        # robust geometric predicates (P1.4a)
katai_add_test(test_delaunay katai::mesh)          # incremental Delaunay (P1.4b)
katai_add_test(test_mohr_coulomb katai::materials)      # MC constitutive kernels (P1.2a)
katai_add_test(test_flux_boundary katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-FLW-002: the prescribed boundary flux (positive into the soil, flow problem only) against Darcy's straight head profile, with the discharge closing through the base
katai_add_test(test_dilatancy_cutoff katai::materials)  # KV-CST-004: the dilatancy cut-off -- dilation stops at the critical void ratio, the flow is isochoric afterwards, and the switch is inert until it acts
katai_add_test(test_hardening_soil katai::materials)    # HS stiffness laws + triaxial hyperbola (P2.3a)
katai_add_test(test_hs_shear katai::materials)          # HS shear return mapping, triaxial (P2.3d-1)
katai_add_test(test_hs_cap katai::materials)            # HS cap return mapping, isotropic (P2.3d-2)
katai_add_test(test_hs_fe katai::materials)             # HS integrate_point FE binding (P2.3d-3a)
katai_add_test(test_hs_two_surface katai::materials)    # HS shear+cap two-surface return (P2.3d-2b)
katai_add_test(test_hs_integrate katai::materials)      # HS robust substepping + Koiter integrator (P2.3d-2b)
katai_add_test(test_hs_berlin katai::materials)         # HS Berlin Sand III drained triaxial validation (contraction then dilation, MC plateau)
katai_add_test(test_hs_undrained_triaxial katai::materials)  # HS Berlin Sand III UNDRAINED triaxial
katai_add_test(test_hs_undrained_clay katai::materials)  # HS undrained NC clay (psi=0), effective MC envelope, Undrained (A) effective-stress analysis
katai_add_test(test_hs_calibration katai::materials)    # HS cap calibration K0nc+Eoed (von Mises, P2.3d-2b)
katai_add_test(test_plate katai::fem)             # Plate (Mindlin/Timoshenko beam) cantilever validation (P2.4)
katai_add_test(test_plate5 katai::fem)            # 5-node quartic plate (tri15 edge) validation (P2.4)
katai_add_test(test_consolidation katai::analysis)     # Biot coupled consolidation -- Terzaghi 1D U-Tv (Phase A.2)
katai_add_test(test_transient_flow katai::analysis katai::mesh)  # Transient saturated flow (W1): Terzaghi isochrones + Carslaw-Jaeger erfc + steady regression
katai_add_test(test_water_retention katai::materials)   # van Genuchten SWCC + Mualem k_rel (W2a): limits/monotonicity/asymptotes/Mualem n=2 + analytic capacity = FD
katai_add_test(test_unsaturated_flow katai::analysis katai::mesh)  # Transient unsaturated Richards (W2b): saturated Terzaghi + Celia mass conservation + Philip sqrt(t)
katai_add_test(test_coupled_flow katai::analysis katai::mesh)  # Fully-coupled flow-deformation (W3): saturated = consolidation bit-for-bit + Terzaghi + unsaturated Bishop activation
katai_add_test(test_dynamics katai::analysis)          # Dynamics core (seismic D1): Newmark-beta + Rayleigh -- exact SDOF + 1D site response f_n=(2n-1)Vs/4H + base-harmonic transfer function
katai_add_test(test_material_profile katai::analysis katai::mesh katai::geometry)  # Depth-varying E(y)/c(y) (v0.5): uniform limit BIT-FOR-BIT + closed form q*integral dy/Eoed(y) + convergence of a constant-E sublayer stack + same profile on the nonlinear path + c_inc reaching strength
katai_add_test(test_ssi_dynamics katai::analysis katai::mesh katai::geometry)  # Soil-structure interaction in a dynamic system (D5): structural K == static solve_nonlinear tangent (machine precision) + K.r=0 rigid mode + r'Mr = total mass + plate cantilever eigenfrequency (Euler-Bernoulli)
katai_add_test(test_dynamics_nonlinear katai::analysis katai::mesh katai::geometry)  # Nonlinear dynamics (Newmark + per-step Newton, v0.5): linear limit -> solve_newmark EXACTLY + quasi-static limit -> static elastoplastic solve_nonlinear (MC yielding) + a witness of nonlinearity
katai_add_test(test_tbdy_seismic katai::analysis)      # TBDY 2018 seismic (D4a): horizontal elastic design spectrum + site coefficients Tables 2.1/2.2 (official AFAD) + response spectrum from an accelerogram (SDOF solve_newmark)
katai_add_test(test_consolidation_plastic katai::analysis)  # Elastoplastic (MC/HS) Biot consolidation: LE reduction = Terzaghi + MC coupled Newton
katai_add_test(test_coupled_flow_plastic katai::analysis katai::mesh)  # Elastoplastic (MC/HS) fully-coupled flow-deformation (W3 follow-up): saturated = consolidation_plastic bit-for-bit + LE unsaturated = W3 + MC unsaturated robustness
katai_add_test(test_general_interface katai::analysis katai::mesh)  # General interface (arbitrary orientation, soil-soil): rigid = bonded continuum + soft slip (present != absent), tri6+tri15
katai_add_test(test_hssmall katai::materials)           # HSsmall small-strain stiffness law (Phase A.5)
katai_add_test(test_material_axisym_coverage katai::materials) # every registered model is integrated in BOTH modes: outputs poisoned with NaN first, so a switch branch that does not exist cannot pass by handing back the caller's own values
katai_add_test(test_material_registry katai::materials) # constitutive catalogue: name resolution, build == the driver's former inline construction (field-exact), honest drainage refusals, external registration through the one seam
katai_add_test(test_design_code katai::analysis)       # EC7 (EN 1997-1) + TBDY 2018 partial factor tables + material factoring (v0.3 B1)
katai_add_test(test_geometry2d katai::geometry)     # planar arrangement: geometry lines -> closed regions (clusters)
katai_add_test(test_model_mesh katai::jobs katai::mesh katai::model)  # project polygon -> FE mesh (katai/jobs/mesh_builder.hpp)
katai_add_test(test_plate_mesh katai::jobs katai::mesh katai::model)  # regression: a structural line does not delete the mesh (outline separation)
katai_add_test(test_mesh_domain katai::jobs katai::mesh katai::model)  # mesh triangles must NOT spill outside the domain (load/structure/non-convex)
katai_add_test(test_mesh_quality_signal katai::jobs katai::mesh katai::model)  # KV-DIA-002: a mesh says whether it reached its own quality bound. Ruppert guarantees the angle, but the loop carries a step cap and the valve used to open in silence -- the K2D-A012 argument, applied to the mesher
katai_add_test(test_k0_layered katai::analysis katai::mesh)  # general layered + water K0 initial stress
katai_add_test(test_point_location katai::fem)    # Point location + inverse isoparametric map (embedded beam, Phase A.4)
katai_add_test(test_embedded_beam katai::analysis)     # Embedded beam skin interaction -- axial pile closed form (Phase A.4)
katai_add_test(test_distributed_load katai::fem katai::mesh)  # Distributed (line) load -- consistent varying-traction assembly
katai_add_test(test_struct_forces katai::analysis katai::mesh)  # Structural internal-force output: plate M/Q/N diagram + anchor/geogrid
katai_add_test(test_assembly katai::fem)   # pipeline: assembly -> solve -> analytic (P0.6)
katai_add_test(test_stress katai::analysis)      # stress recovery (P0.7+P0.8a)
katai_add_test(test_boussinesq katai::analysis)  # Boussinesq strip load (P0.8b)
katai_add_test(test_flamant katai::analysis)     # Flamant point/line load sigma_z=2P/pi.z (Phase A.3)
katai_add_test(test_nonlinear katai::analysis)   # Newton-Raphson skeleton (P1.1)
katai_add_test(test_prandtl katai::analysis)     # Prandtl Nc=5.14 limit load (phi=0, P1.3)
katai_add_test(test_bearing_phi katai::analysis) # Reissner Nc(phi=20deg) bearing capacity (P1.3)
katai_add_test(test_mesher_pipeline katai::analysis)  # mesher -> tri6 -> FEM (P1.4e)
katai_add_test(test_tri15_pipeline katai::fem)   # mesher -> tri15 -> FEM end to end (P1.5)
katai_add_test(test_tri15_traction katai::fem)   # tri15 quartic-edge surface traction (P1.5)
katai_add_test(test_tri15_prandtl katai::analysis)    # tri15 Prandtl Nc (coarse mesh -> 1.1%, P1.5)
katai_add_test(test_axisym_cylinder katai::fem)  # axisymmetric Lame cylinder (P1.6)
katai_add_test(test_axisym_collapse katai::analysis)  # axisymmetric MC cylinder collapse 2c.ln(b/a) (P1.6)
katai_add_test(test_effective_stress katai::fem) # effective stress + pore pressure (water)
katai_add_test(test_phreatic_gravity katai::fem) # phreatic (water table) buoyant gravity + effective stress
katai_add_test(test_undrained katai::materials)        # undrained (A) constitutive + Skempton B
katai_add_test(test_undrained_global katai::analysis) # global undrained (A) confined loading
katai_add_test(test_initial_stress_k0 katai::analysis) # K0 procedure + gravity loading (P2.2)
katai_add_test(test_excavation katai::analysis)        # staged construction excavation heave (P2.2)
katai_add_test(test_staged_fill katai::analysis)       # staged fill + multi-phase chain (P2.2)
katai_add_test(test_hs_fe_oedometer katai::analysis)   # HS two-surface (cap+shear) FE/BVP validation (P2.3)
katai_add_test(test_hs_footing katai::analysis katai::mesh)  # HS strip footing settlement (shear-dominated BVP robustness)
katai_add_test(test_hs_undrained_bvp katai::analysis katai::mesh)  # HS + Undrained(A) BVP (clay footing, wrapper+HS)
katai_add_test(test_hs_axisym katai::analysis katai::mesh)   # HS axisymmetric (oedometer K0 + circular footing BVP)
katai_add_test(test_plate_soil katai::analysis)        # Plate-in-soil coupling (DofMap rotation + solver, P2.4)
katai_add_test(test_plate_plastic katai::analysis)     # Plate M-N hinge (Mp/Np diamond interaction, perfectly plastic)
katai_add_test(test_anchor katai::analysis)            # Anchor (fixed-end + node-to-node axial spring, P2.4)
katai_add_test(test_anchor_plastic katai::analysis)    # Anchor elastoplastic F_max + permanent U_p (P2.4)
katai_add_test(test_geogrid katai::analysis)           # Geogrid tension-only axial membrane + N_p (P2.4)
katai_add_test(test_interface katai::analysis)         # Interface (Goodman zero-thickness Coulomb, P2.4)
katai_add_test(test_seepage katai::analysis)           # Steady-state confined seepage (Darcy/Laplace, P3)
katai_add_test(test_earth_pressure katai::analysis)    # Rankine active/passive lateral earth pressure (P2.4 excavation)
katai_add_test(test_wall_split katai::analysis)        # Embedded wall mesh node splitting (barrier, P2.4)
katai_add_test(test_embedded_wall katai::analysis)    # embedded wall cantilever: mesh split + plate extra DOFs + two interfaces
katai_add_test(test_excavation_wall katai::analysis)  # Coupled embedded wall + soil excavation (bonded plate-in-soil, P2.4)
katai_add_test(test_wall_k0_excavation katai::analysis) # Barrier wall + K0 interface seeding + excavation (P2.4)
katai_add_test(test_slope katai::analysis)            # phi-c reduction slope FoS (P1.7)
katai_add_test(test_solve katai::jobs katai::mesh katai::analysis katai::model)  # project -> mesh -> LE self-weight solve (katai/jobs/driver.hpp)
katai_add_test(test_gui_solve katai::jobs katai::mesh katai::analysis katai::model)  # GUI solve K0 + water + effective stress, end to end
katai_add_test(test_anchor_mesh_repro katai::jobs katai::mesh katai::analysis katai::model)  # regression: an anchor does not corrupt the mesh (nearest node)
katai_add_test(test_constitutive katai::jobs katai::mesh katai::analysis katai::model)  # the GUI constitutive model is genuinely in effect, LE/MC/HS
katai_add_test(test_hoek_brown katai::materials)   # Hoek-Brown (rock): the closed forms of Hoek, Carranza-Torres & Corkum (2002), the envelope, and the Tresca degeneration against the verified MC return (KV-CST-014)
katai_add_test(test_hoek_brown_tunnel katai::jobs katai::mesh katai::analysis katai::model)  # a tunnel unloaded into Hoek-Brown rock, against the closed form (KV-CST-016): the case that found the missing axisymmetric branch, the elastic tangent and the retry that never reached it
katai_add_test(test_hoek_brown_gui katai::jobs katai::mesh katai::analysis katai::model)  # Hoek-Brown from a .k2d: a confined specimen fails ON the generalised Hoek-Brown envelope, at two confinements (KV-CST-015)
katai_add_test(test_axisym_gui katai::jobs katai::mesh katai::analysis katai::model)  # Axisymmetric GUI path: K0 + r-weighted gravity/traction + radial hoop
katai_add_test(test_axisym_water katai::jobs katai::mesh katai::analysis katai::model katai::fem katai::geometry)  # Water table in axisymmetric mode: the K0 equilibrium identity and the HOOP pore term (KV-CST-013)
katai_add_test(test_gui_matrix katai::jobs katai::mesh katai::analysis katai::model)  # Systematic: structures x load x water x model x tri6/tri15 invariants
katai_add_test(test_safety_gui katai::jobs katai::mesh katai::analysis katai::model)  # Safety (phi-c reduction) GUI path: slope FoS + collapse mechanism
katai_add_test(test_k0_slope katai::jobs katai::mesh katai::analysis katai::model)    # K0 nil-step equilibrium on a sloping surface (horizontal = exact u0; a real redistribution on the slope)
katai_add_test(test_seepage_gui katai::jobs katai::mesh katai::analysis katai::model) # Groundwater flow GUI path: exact Darcy + Charny seepage face + Terzaghi coupling
katai_add_test(test_consolidation_gui katai::jobs katai::mesh katai::analysis katai::model) # Consolidation GUI path (PhaseType::Consolidation): 1D Terzaghi U-Tv + final settlement
katai_add_test(test_consolidation_stop katai::jobs katai::mesh katai::analysis katai::model katai::linsolve) # Consolidation stop criteria (phases[i].cstop, docs/k2d-format.md): restart identity + the PRESSURE-ratio degree of consolidation against Terzaghi
katai_add_test(test_consolidation_structures katai::jobs katai::mesh katai::analysis katai::model) # Plates/anchors/geogrids INSIDE a consolidation phase: the drained limit is the oracle (KV-STR-006)
katai_add_test(test_consolidation_layers katai::jobs katai::mesh katai::analysis katai::model)  # KV-CON-004: the storage term n/Kw per material -- each layer's undrained pressure, the same answer whatever order the materials are listed in, and the settlement history against an independent finite-volume column
katai_add_test(test_consolidation_interface katai::jobs katai::mesh katai::analysis katai::model)  # Interfaces in a consolidation phase: cross permeability ties/separates the pore field, elastic joint reported (KV-STR-007)
katai_add_test(test_fully_coupled_structures katai::jobs katai::mesh katai::analysis katai::model)  # Structures in a FULLY-COUPLED phase: the saturated limit IS the consolidation phase (KV-STR-008)
katai_add_test(test_dynamic_gui katai::jobs katai::mesh katai::analysis katai::model) # Dynamic/seismic GUI path (PhaseType::Dynamic, D4b): 1D site-response resonance (4/pi A/w1^2 2xi) + frequency selectivity
katai_add_test(test_water_gui katai::jobs katai::mesh katai::analysis katai::model) # Water-history GUI path: FullyCoupled = Terzaghi + TransientFlow flow (W4 wiring)
katai_add_test(test_undrained_gui katai::jobs katai::mesh katai::analysis katai::model) # Undrained (A/B) GUI path: confined column -qH/M_u + effective stress (Skempton 1D)
katai_add_test(test_robustness katai::jobs katai::mesh katai::analysis katai::model) # Robustness: the calculation path never crashes or leaks NaN on hostile/degenerate input
katai_add_test(test_struct_gui katai::jobs katai::mesh katai::analysis katai::model)  # Structural M/Q/N diagrams from the GUI path (wiring: stations along the line, free-end M0, strut N)
katai_add_test(test_interface_gui katai::jobs katai::mesh katai::analysis katai::model)  # Standalone interface GUI path: embed in mesh + split + soil-soil Coulomb -> present != absent (tri6+tri15)
katai_add_test(test_mesh_refine katai::jobs katai::mesh katai::model) # Local mesh density (coarseness/sizing field): region + line refinement, quality, bit-for-bit legacy
katai_add_test(test_staged_gui katai::jobs katai::mesh katai::analysis katai::model)  # Multi-phase staged construction (each phase ramps its configuration imbalance): fill settlement + excavation heave + phase-chain load (closed form)
katai_add_test(test_staged_struct_carry katai::jobs katai::mesh katai::analysis katai::model)  # Static chain structural state carry: NIL phase no-op IDENTITY (round-off; caught a 32% M shift on 2026-07-19 -> fixed)
katai_add_test(test_site_response_benchmark katai::jobs katai::mesh katai::analysis katai::model)  # BENCHMARK Wave-1: multi-layer site response vs a layered viscoelastic transfer-matrix solution (independent oracle, Track 7)
katai_add_test(test_compliant_base katai::jobs katai::mesh katai::analysis katai::model)  # Compliant (absorbing) base vs a layered solution with a radiation BC: uniform |T|~1, rigid-vs-compliant tooth, honest guards
katai_add_test(test_real_record katai::jobs katai::mesh katai::analysis katai::model)  # Real record: accelerogram input (bit-for-bit wiring) + El Centro 1940 NS (published PGA + spectral shape) + a product run
katai_add_test(test_foundation_benchmarks katai::mesh katai::analysis katai::geometry)  # KV-FND-001..004: footing benchmarks against analytical solutions -- Giroud (1972) rigid footing, Gibson (1967) E(z), Cox (1962) circular bearing, Davis & Booker (1973) c(z)
katai_add_test(test_soft_soil katai::mesh katai::analysis katai::geometry)  # Soft Soil: stage 1 material point (ln law lambda*/kappa*, memory, K0NC pin, MC failure) + stage 2 FE wiring (integrate_point identity, frame independence, axisymmetry, oedometer BVP)
katai_add_test(test_soft_soil_gui katai::jobs katai::mesh katai::analysis katai::model)  # Soft Soil stage 3 GUI path: K0 + ss_initial_pp seeding, NC lambda* integral, OCR/POP memory, honest undrained refusal, dynamic p'/kappa* stiffness
katai_add_test(test_soft_soil_creep katai::mesh katai::analysis katai::geometry)  # Soft Soil Creep: stage 1 material point (mu*.ln(1+t/tau), OCR^-beta, closed-form relaxation, 1-day NC line, MC) + stage 2 FE (time into integrate_point, time-stepped relaxation BVP identity) (Vermeer & Neher 1999)
katai_add_test(test_input_audit katai::jobs katai::mesh katai::analysis katai::model)  # Input audit pins: OCR -> K0 = K0nc*OCR - nu/(1-nu)*(OCR-1), HS+Undrained(B) refusal, plate Mp/Np refusal
if(TARGET test_real_record)
  target_compile_definitions(test_real_record PRIVATE KATAI_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
endif()
katai_add_test(test_design_gui katai::jobs katai::mesh katai::analysis katai::model)  # EC7 design approach through the GUI path (build_problem): DA3 ODF = FoS/gamma_M identity + None baseline (v0.3 B1.3)
katai_add_test(test_project_io katai::model katai::io)                # Project file JSON round trip (every field exact + honest failure modes)
katai_add_test(test_results_io katai::mesh katai::analysis katai::model katai::io)  # Binary result file: bitwise round trip + model-hash staleness + refusal of a corrupt file
katai_add_test(test_k2d_validator katai::model katai::io)             # .k2d input contract: every rule probed at its exact field path + reader forward-version notes
katai_add_test(test_stratigraphy katai::model)                        # KV-GEO-001: borehole logs -> soil regions + water surface (strata/boreholes, docs/k2d-format.md). A wrong stratigraphy is invisible downstream -- the mesher meshes it, every phase converges on it, and no residual sees it
# One ctest entry per corpus case. Together they were 1370 s and, once test_phase_numerics had
# been split, the whole suite's wall clock; they are independent by construction, so this is the
# same trade the phase-numerics split made. The ids are the corpus ids, so a red line in the ctest
# output names the case without anyone opening the log. The executable still runs all 26 when
# invoked with no argument -- which is also what KATAI_CORPUS_WRITE=1 regeneration needs.
set(KATAI_TEST_CASES
  kv-con-002 kv-cst-001 kv-cst-002 kv-cst-008 kv-cst-009 kv-cst-010 kv-cst-011 kv-cst-012
  kv-dyn-002 kv-dyn-003 kv-exc-001 kv-flw-001
  kv-fnd-008 kv-fnd-009 kv-fnd-010 kv-fnd-011 kv-fnd-012 kv-fnd-013 kv-fnd-014
  kv-num-003 kv-slp-001 kv-slp-002 kv-str-002 kv-str-003 kv-str-004 kv-str-005)
list(LENGTH KATAI_TEST_CASES _katai_corpus_ids)
katai_add_test(test_input_corpus katai::jobs katai::mesh katai::analysis katai::model katai::io)  # Input corpus (roadmap 6.3): checked-in .k2d byte-identical to the programmatic build, validated, then SOLVED FROM THE FILE against its declared closed form
unset(KATAI_TEST_CASES)
if(TARGET test_input_corpus)
  target_compile_definitions(test_input_corpus PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus"
    KATAI_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"    # El Centro record (KV-DYN-003)
    # ...and the count the ids above registered, which the test static_asserts against its own
    # table. See the comment there: without it a case can exist with no ctest entry.
    KATAI_CORPUS_CASE_COUNT=${_katai_corpus_ids})
endif()
katai_add_test(test_phase_water katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-CST-003: staged dewatering (.k2d v4 per-phase water) against the effective-stress closed form -- the second capability the anchored-excavation benchmark needed
katai_add_test(test_anchor_prestress katai::jobs katai::mesh katai::analysis katai::model)  # KV-STR-001: the anchor lock-off force (.k2d v3) against an equivalent external-force model -- the first structural case in the matrix, and the input every anchored excavation needs
katai_add_test(test_driven_structures katai::jobs katai::mesh katai::analysis katai::model)  # KV-STR-009: a plate on a settling support and anchors on driven edges, read from the force report against closed forms -- the anchor report used to drop a driven node's motion
katai_add_test(test_safety_structures katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-STR-010: a Safety run with its structures, on the sliding block -- the interface reduced with the soil, a plate's weight and an anchor's capacity carried unreduced, each against a closed form
if(TARGET test_safety_structures)
  target_compile_definitions(test_safety_structures PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
katai_add_test(test_tolerance_independence katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-NUM-007: the same benchmark at six tolerated residuals -- is the published number a physics result or an artefact of the stopping rule? (and, above the default, an unsafe-sided one)
if(TARGET test_tolerance_independence)
  target_compile_definitions(test_tolerance_independence PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
katai_add_test(test_convergence_criteria katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-NUM-013: the convergence criteria FAMILY -- CSP, the CSP-normalised force error and the local error at every soil stress point; the elastic control that must report nothing, and what it costs to be bound by criteria this tree measured but did not obey
if(TARGET test_convergence_criteria)
  target_compile_definitions(test_convergence_criteria PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
katai_add_test(test_flow_barrier katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-FLW-004: a wall as a flow barrier -- an impermeable screen passes exactly no water and holds the whole head difference, a semi-permeable one passes dh / (L/(k b) + R/b)
katai_add_test(test_wells_drains katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-FLW-003: wells and drains against the Dupuit 1D closed forms -- the drawdown of a full-depth well, the discharge a drain removes, and a well limited by h_min becoming a drain at that level
katai_add_test(test_undrained_c katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-CST-007: Undrained (C), the total-stress analysis -- the same soil described two ways must agree, and the difference is the compressibility nu_u < 0.5 leaves behind
katai_add_test(test_undrained_stiffness katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-CST-006: the pore fluid's stiffness per material (Kw/n, nu_u or Skempton's B) -- two undrained layers with different Kw/n, and the Hardening Soil family sized by its own unload/reload pair
katai_add_test(test_staged_options katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-EXC-002 + KV-CST-005: the staged-construction fraction (mstage) as an exact linear identity, and "ignore undrained behaviour" as the drained calculation it claims to be
# FOUR independent verification cases in one file, so four ctest entries: KV-NUM-008 (the
# per-phase numerical controls of .k2d v7 must REACH the solver -- file == seam, and both !=
# default, so a reviewer receives the numerics a published run was computed with), KV-NUM-009
# (where this case's numerical error actually lives), KV-NUM-010 (the same procedure on an axis
# that IS clean) and KV-NUM-011 (where "the algebra decides the order" stops being true). Run
# together they were 2683 s -- measured 2026-08-26 -- which was the entire suite's wall clock,
# because `ctest -j 6` cannot finish faster than its longest single test.
set(KATAI_TEST_CASES 008 009 010 011)
katai_add_test(test_phase_numerics katai::jobs katai::mesh katai::analysis katai::model katai::io)
unset(KATAI_TEST_CASES)
if(TARGET test_phase_numerics)
  target_compile_definitions(test_phase_numerics PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
# Hardening Soil stress-point integration: the instrument KV-NUM-009 was re-measured with
# (EXCLUDE_FROM_ALL, like the other studies; deliberately not add_test'd).
if(TARGET katai::jobs AND TARGET katai::io)
  add_executable(study_hs_integration EXCLUDE_FROM_ALL study_hs_integration.cpp)
  target_link_libraries(study_hs_integration PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::io katai::linsolve)
  target_compile_definitions(study_hs_integration PRIVATE KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")

  # The convergence criteria family (N-2): prints every member for a set of cases, so the
  # question "does a case exist that meets the global criterion and fails a local one?" is
  # answered by a run rather than by an argument.
  add_executable(study_convergence_family EXCLUDE_FROM_ALL study_convergence_family.cpp)
  target_link_libraries(study_convergence_family PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::io katai::linsolve)
  target_compile_definitions(study_convergence_family PRIVATE KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()

katai_add_test(test_slope_fos_convergence katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-SLP-003: the strength-reduction factor of safety of the Griffiths-Lane benchmark under systematic nested refinement -- does it converge, or does it follow the shear band down the mesh?
if(TARGET test_slope_fos_convergence)
  target_compile_definitions(test_slope_fos_convergence PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
katai_add_test(test_mesh_refine_uniform katai::jobs katai::mesh katai::model)  # KV-NUM-006: regular ("red") refinement -- four similar children per triangle, nested nodes, angles and area preserved exactly, so a convergence study varies the mesh SIZE and nothing else
katai_add_test(test_mesh_convergence katai::jobs katai::mesh katai::analysis katai::model katai::io katai::fem)  # KV-NUM-005: the strip-load benchmark at three systematically refined densities -> observed order, Richardson extrapolation and the GCI band that every published number is meant to carry
if(TARGET test_mesh_convergence)
  target_compile_definitions(test_mesh_convergence PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
# One ctest entry per fixture: fourteen models built and solved, sharing nothing but the failure
# counter, and 1071 s of wall clock that `ctest -j 6` could not overlap with anything else.
set(KATAI_TEST_CASES
  reference load_above_surface load_outside_model load_half_outside
  point_load_off_mesh point_load_between_nodes plate_above_soil plate_clipped
  geogrid_above_soil anchor_outside_soil wall_above_soil tension_cutoff_ignored
  prescribed_disp_on_structure linear_dynamic_zero_stress)
list(LENGTH KATAI_TEST_CASES _katai_drop_ids)
katai_add_test(test_no_silent_drop katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-DIA-001: no input may be discarded in silence -- one fixture per driver build loop, each perturbing ONE field of a corpus case, asserting the stable diagnostic code and that a clipped object equals the object drawn where the mesh took it
unset(KATAI_TEST_CASES)
if(TARGET test_no_silent_drop)
  target_compile_definitions(test_no_silent_drop PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus"
    KATAI_DROP_CASE_COUNT=${_katai_drop_ids})
endif()
katai_add_test(test_jobs katai::jobs katai::mesh katai::analysis katai::model katai::io)  # The Job object (Stage D): driver-path equivalence bit-identical, progress, cooperative phase-boundary cancel, timings, two concurrent jobs in one process
if(TARGET test_jobs)
  target_compile_definitions(test_jobs PRIVATE KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
endif()
katai_add_test(test_api_surface katai::api katai::mesh katai::analysis katai::model katai::io)  # The published facade (Stage D): one header suffices for construct/save/load/validate/run/read, and the job refuses an invalid project at its field path
katai_add_test(test_cli katai::api katai::mesh katai::analysis katai::model katai::io)  # The katai command line (7.1): documented exit codes, the corpus as its test set, refusal at field path, CROSS-PROCESS bit-identical .res vs an in-process Job
if(TARGET test_cli AND TARGET katai_cli)
  add_dependencies(test_cli katai_cli)
  target_compile_definitions(test_cli PRIVATE
    KATAI_BIN_DIR="${CMAKE_BINARY_DIR}/bin"
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
elseif(TARGET test_cli)
  # Without the CLI target the test cannot run; exclude it rather than fake it.
  set_tests_properties(test_cli PROPERTIES DISABLED TRUE)
endif()

# PARDISO's own adapter tests: these name the backend on purpose, so they are the
# one thing that legitimately requires MKL. Everything else above runs on whichever
# backend the configuration selected.
if(KATAI_WITH_MKL AND MKL_FOUND)
  katai_add_test(test_pardiso_smoke)        # raw MKL API link/correctness smoke
  katai_add_test(test_pardiso_solver)       # PardisoSolver wrapper unit test
endif()

# --- Validation studies (not regressions) ----------------------------------
# EXCLUDE_FROM_ALL, so they never slow the normal build, and not registered with
# CTest. Run one with `cmake --build <dir> --target <name>`. Results are written
# up under docs/validation/. The ones that drive the GUI compute path need the
# project schema, so they exist only in compositions that build katai::model.
add_executable(study_giroud_convergence EXCLUDE_FROM_ALL study_giroud_convergence.cpp)  # solution verification of the Giroud benchmark: nested refinement + GCI band -> the published record
target_link_libraries(study_giroud_convergence PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::io katai::linsolve)
add_executable(study_slope_convergence EXCLUDE_FROM_ALL study_slope_convergence.cpp)  # mesh convergence of the slope SRM -> docs/validation/slope-srm-convergence.md
target_link_libraries(study_slope_convergence PRIVATE katai::analysis katai::linsolve)
add_executable(study_design_ec7 EXCLUDE_FROM_ALL study_design_ec7.cpp)  # EC7 DA3 material factoring (ODF=FoS/gamma_M identity + Griffiths & Lane benchmark)
target_link_libraries(study_design_ec7 PRIVATE katai::analysis katai::mesh katai::linsolve)
add_executable(study_wall_benchmark EXCLUDE_FROM_ALL study_wall_benchmark.cpp)  # cantilever sheet pile in undrained clay: domain-size audit of the wall moment
target_link_libraries(study_wall_benchmark PRIVATE katai::analysis katai::linsolve)
add_executable(study_circular_footing EXCLUDE_FROM_ALL study_circular_footing.cpp)  # rigid circular footing on sand (axisymmetric MC + K0 + water table), force and performance
target_link_libraries(study_circular_footing PRIVATE katai::analysis katai::linsolve)
add_executable(study_hs_calibration EXCLUDE_FROM_ALL study_hs_calibration.cpp)  # HS alpha, beta calibration
target_link_libraries(study_hs_calibration PRIVATE katai::materials katai::linsolve)
# LBB / inf-sup diagnostic (audit P2b): does equal-order u-p Biot checkerboard under 2D undrained?
add_executable(study_lbb_undrained EXCLUDE_FROM_ALL study_lbb_undrained.cpp)
target_link_libraries(study_lbb_undrained PRIVATE katai::mesh katai::analysis katai::geometry katai::linsolve)
if(TARGET katai::model)
  add_executable(study_default_convergence EXCLUDE_FROM_ALL study_default_convergence.cpp)  # v0.4 diagnosis: is the collapse weak-soil physics? (self-weight converges, the footing collapses at weak c)
  target_link_libraries(study_default_convergence PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::linsolve)
  add_executable(diag_k0 EXCLUDE_FROM_ALL diag_k0.cpp)  # K0 procedure accuracy audit (GUI compute path)
  target_link_libraries(diag_k0 PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::linsolve)
  # Quantitative V&V of the GUI path (on demand, ~45 s): MC Prandtl/Reissner bearing
  # capacity + the elastic regime. The nonlinear collapse solve is expensive, so it is
  # a study rather than a test. Results: docs/validation/gui-pipeline-validation.md.
  add_executable(study_gui_validation EXCLUDE_FROM_ALL study_gui_validation.cpp)
  target_link_libraries(study_gui_validation PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::linsolve)
  # Performance baseline (on demand): where the time goes on the GUI compute path
  # (assembly/CSR/solve). Results: docs/validation/performance-baseline.md.
  add_executable(study_perf EXCLUDE_FROM_ALL study_perf.cpp)
  target_link_libraries(study_perf PRIVATE katai::jobs katai::mesh katai::analysis katai::model katai::linsolve)
endif()

# --- Source-tree gates (Python; no compilation, so they cost seconds) -------
# check_language: KATAI is meant to be read and audited internationally, so every
#   comment and identifier in the tree must be English. The translation sweep
#   spans many commits, so the gate runs in ratchet mode - existing debt is
#   tolerated, but a file that GAINS non-English prose fails the build.
# check_architecture: the layer contract (engine -> jobs -> api -> front ends)
#   only stays true if something enforces it. The gate reads the include graph
#   and refuses upward dependencies, cycles, private-header reach-through and
#   physics defined above the engine. Known violations are listed in the script
#   with the stage that removes them, so an exception is a debt with an owner.
# Both are source analysis rather than C++ tests on purpose: adding translation
# units to a tree whose clean build already costs an hour would be self-defeating.
# One FindPython module across the whole tree: nanobind (python/) requires the NEW
# module, and mixing FindPython with FindPython3 corrupts each other's cached state
# (measured: a reconfigure after adding the bindings silently dropped every gate).
find_package(Python COMPONENTS Interpreter QUIET)
if(Python_Interpreter_FOUND)
  add_test(NAME check_language
           COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/check_language.py --ratchet)
  set_tests_properties(check_language PROPERTIES LABELS "hygiene")
  add_test(NAME check_architecture
           COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/check_architecture.py)
  set_tests_properties(check_architecture PROPERTIES LABELS "architecture")
  # test_reference_registry (roadmap section 6.2): every verification test declares
  # its oracle -- source, locator, expected value and band -- or is in the frozen
  # adoption list; the verification matrix and references.bib must equal their
  # regeneration from those declarations. Named per the roadmap; implemented as a
  # source gate because it reads comments, not binaries.
  add_test(NAME test_reference_registry
           COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/check_references.py)
  set_tests_properties(test_reference_registry PROPERTIES LABELS "references")
  # check_k2d_spec (Stage C): the .k2d format's four descriptions -- writer,
  # reader, docs/k2d-format.md field inventory, docs/k2d.schema.json -- must
  # agree on the key set, the version numbers and the file-stable enum bounds.
  # A field added to the code without its documentation fails here, not in a
  # user's hands.
  add_test(NAME check_k2d_spec
           COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/check_k2d_spec.py)
  set_tests_properties(check_k2d_spec PROPERTIES LABELS "io")
  # check_silent_drop: the safety property test_no_silent_drop asserts at runtime, held
  # at the source. In the driver's build loops -- where a drawn object becomes an engine
  # object -- a `continue` that skips user input must emit a diagnostic, be a selection
  # (wrong kind, not in this phase), or carry a stated reason. A test can only cover the
  # cases someone thought of; this covers the ones nobody has written yet.
  add_test(NAME check_silent_drop
           COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/check_silent_drop.py)
  set_tests_properties(check_silent_drop PROPERTIES LABELS "hygiene")
else()
  message(STATUS "Python interpreter not found - source-tree gates (language, architecture) disabled")
endif()

# THE SLOW SET, measured rather than guessed, so that the inner development loop does not
# have to pay for the full suite. Measured 2026-08-24 on the MKL composition (ctest -j 6):
# the 154 tests cost 9089 s of CPU between them and TWELVE of them are 8461 s of that -- 93%.
# With `-j 6` the wall clock cannot fall below the single longest test, so the suite takes as
# long as its slowest file does.
#
#   ctest -LE slow -j 6     the inner loop: everything except these, ~600 s of CPU
#   ctest -j 6              before a commit: everything
#
# The second number on each line is the measured cost in seconds, given to CTest so that it
# schedules the long ones FIRST on a cold checkout too (CTest learns costs from previous runs,
# but only after there has been one). A number here going badly stale is visible in the ctest
# output as poor packing, not as a wrong answer, which is the right failure mode for a hint.
#
# 2026-08-26, three splits in one day, each of a file that was the suite's wall clock at the time:
#
#   test_phase_numerics   2683 s -> four entries, worst 1176 s   (four independent KV cases)
#   test_input_corpus     1370 s -> 26 entries, worst  306 s     (one per corpus file)
#   test_no_silent_drop   1071 s -> 14 entries, worst  469 s     (one per dropped-input fixture)
#
# AND THAT IS AS FAR AS SPLITTING GOES. Splitting buys wall clock only while one entry exceeds
# the suite's total divided by N, and none does any more: the longest is test_hs_footing at 890 s
# solo against a 1796 s run. Past that the lever is the work the solves do, not the size of the
# entries -- which is exactly what the next line shows, since 2026-08-27's 470 s came from ONE
# verification block giving up a rung it should never have been standing on (KV-NUM-009 (d)):
#
#   full suite, MKL/PARDISO      154 tests / 2597 s  ->  195 / 2266 s  ->  195 / 1796 s
#   full suite, portable/Eigen                           193 /  ---   ->  193 / 1881 s
#   inner loop                   144 tests /  144 s  ->  180 /  169 s
#
# The gain beyond the wall clock is that 36 real verification cases fell below the slow threshold
# and joined the INNER loop, where before they were checked only in the pre-commit run.
#
# OVER-SUBSCRIPTION IS REAL AND COSTS NOTHING, WHICH IS WHY `-j 6` STAYS. Six slots on four
# physical cores does slow each heavy solve down, and the two full runs say by how much -- these
# are each one run's own accumulated test time for the same fifteen entries, out of that run's
# label summary rather than out of any average:
#
#   ctest -j 6    slow set 8491 s of test time    wall 2266 s
#   ctest -j 4    slow set 6546 s of test time    wall 2257 s
#
# Thirty per cent more test time at -j 6, and the same wall clock to nine seconds. What -j 4 buys
# back per test it gives up again in slots, so the machine is saturated on memory bandwidth either
# way and the scheduler is not where anything is left. Effective parallelism is 3.6x on four
# cores; the lever that remains is the work the solves do.
#
# THE SECONDS BELOW ARE SOLO WALL CLOCK, one entry run by itself, and that is a choice worth
# stating. CTest keeps its own estimate in CTestCostData, but that file holds a RUNNING AVERAGE
# over every run it has seen -- the format is `name count cost`, and test_phase_numerics_009 sat
# at 1735.9 s over 5 runs on a day when the entry itself measured 781 s. Reading a number out of
# that file is not reading a measurement. What is here is what a cold tree will see, and CTest
# replaces it with its own average as soon as it has one.
set(_katai_slow_tests
  test_hs_footing 890             test_phase_numerics_009 781
  test_no_silent_drop_tension_cutoff_ignored 469  test_slope_fos_convergence 535
  test_tolerance_independence 409 test_phase_numerics_008 401
  test_hs_axisym 317              test_input_corpus_kv-cst-012 306
  test_foundation_benchmarks 294  test_convergence_criteria 230
  test_gui_matrix 173             test_dynamic_gui 128    test_bearing_phi 125
  test_input_corpus_kv-str-004 90 test_input_corpus_kv-fnd-013 72)
while(_katai_slow_tests)
  list(POP_FRONT _katai_slow_tests _t _cost)
  if(TEST ${_t})
    set_property(TEST ${_t} APPEND PROPERTY LABELS slow)
    set_property(TEST ${_t} PROPERTY COST ${_cost})
  endif()
endwhile()

# Say out loud what this composition left out, so a shrunken test set is a stated
# fact rather than a surprise in the ctest total.
get_property(_katai_skipped GLOBAL PROPERTY KATAI_SKIPPED_TESTS)
if(_katai_skipped)
  list(LENGTH _katai_skipped _n)
  message(STATUS "katai tests: ${_n} skipped (their modules are not in this composition): ${_katai_skipped}")
endif()
# The same for a narrowed suite. A filter is a build-time saving traded against coverage, and
# which side of that trade a tree is on has to be readable in its configure output, not inferred
# from a test count.
get_property(_katai_filtered GLOBAL PROPERTY KATAI_FILTERED_TESTS)
if(_katai_filtered)
  list(LENGTH _katai_filtered _n)
  message(STATUS "katai tests: ${_n} NOT BUILT here -- KATAI_CORE_TESTS_MATCHING="
                 "'${KATAI_CORE_TESTS_MATCHING}' selects a subset of the core suite. "
                 "The full suite lives in the core's own build tree.")
endif()
