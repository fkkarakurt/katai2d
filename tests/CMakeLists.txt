# tests -- unit and regression. No framework: CTest plus assertions is enough, and
# every test is a plain executable that returns non-zero on failure.
#
# Adding a test: katai_add_test(<name> [link targets...]). The default link is
# katai::math. katai::linsolve is added to every test unconditionally, because a
# test executable is a composition root: it is a whole program, and a program that
# solves has to contain a backend. Which backend that is comes from the
# configuration, never from the test -- see kernel/linsolve/CMakeLists.txt.
#
# The named modules do two more jobs. They label the test, so `ctest -L materials`
# runs exactly the tests whose subject is that module. And they scope it: in a
# composition that does not build one of the named modules the test is skipped,
# never linked against a substitute -- a missing module is refused, not imitated
# (composability rule 3).
function(katai_add_test name)
  set(libs ${ARGN})
  if(NOT libs)
    set(libs katai::math)
  endif()
  set(labels "")
  foreach(lib ${libs})
    if(lib MATCHES "^katai::(.+)$")
      if(NOT TARGET ${lib})
        set_property(GLOBAL APPEND PROPERTY KATAI_SKIPPED_TESTS ${name})
        return()
      endif()
      list(APPEND labels ${CMAKE_MATCH_1})
    endif()
  endforeach()
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ${libs} katai::linsolve)
  katai_use_pch(${name})
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS "${labels}")
  # With MKL linked dynamically its runtime has to be locatable at load time. Putting
  # it in the test's own environment rather than expecting the caller's PATH means
  # `ctest` works from any shell, including a bare one and CI.
  if(KATAI_MKL_RUNTIME_DIR)
    set_tests_properties(${name} PROPERTIES
      ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${KATAI_MKL_RUNTIME_DIR}")
  endif()
endfunction()

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
katai_add_test(test_hardening_soil katai::materials)    # HS stiffness laws + triaxial hyperbola (P2.3a)
katai_add_test(test_hs_shear katai::materials)          # HS shear return mapping, triaxial (P2.3d-1)
katai_add_test(test_hs_cap katai::materials)            # HS cap return mapping, isotropic (P2.3d-2)
katai_add_test(test_hs_fe katai::materials)             # HS integrate_point FE binding (P2.3d-3a)
katai_add_test(test_hs_two_surface katai::materials)    # HS shear+cap two-surface return (P2.3d-2b)
katai_add_test(test_hs_integrate katai::materials)      # HS robust substepping + Koiter integrator (P2.3d-2b)
katai_add_test(test_hs_berlin katai::materials)         # HS Berlin Sand III triaxial validation (PLAXIS Fig 15.4)
katai_add_test(test_hs_undrained_triaxial katai::materials)  # HS Berlin Sand III UNDRAINED triaxial (PLAXIS Fig 15.5)
katai_add_test(test_hs_undrained_clay katai::materials)  # HS undrained NC clay (psi=0), effective MC envelope (PLAXIS Undrained-A)
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
katai_add_test(test_material_registry katai::materials) # constitutive catalogue: name resolution, build == the driver's former inline construction (field-exact), honest drainage refusals, external registration through the one seam
katai_add_test(test_design_code katai::analysis)       # EC7 (EN 1997-1) + TBDY 2018 partial factor tables + material factoring (v0.3 B1)
katai_add_test(test_geometry2d katai::geometry)     # planar arrangement: geometry lines -> closed regions (clusters)
katai_add_test(test_model_mesh katai::jobs katai::mesh katai::model)  # project polygon -> FE mesh (katai/jobs/mesh_builder.hpp)
katai_add_test(test_plate_mesh katai::jobs katai::mesh katai::model)  # regression: a structural line does not delete the mesh (outline separation)
katai_add_test(test_mesh_domain katai::jobs katai::mesh katai::model)  # mesh triangles must NOT spill outside the domain (load/structure/non-convex)
katai_add_test(test_k0_layered katai::analysis katai::mesh)  # general layered + water K0 initial stress
katai_add_test(test_point_location katai::fem)    # Point location + inverse isoparametric map (embedded beam, Phase A.4)
katai_add_test(test_embedded_beam katai::analysis)     # Embedded beam skin interaction -- axial pile closed form (Phase A.4)
katai_add_test(test_distributed_load katai::fem katai::mesh)  # Distributed (line) load -- consistent varying-traction assembly
katai_add_test(test_struct_forces katai::analysis katai::mesh)  # Structural internal-force output: plate M/Q/N diagram + anchor/geogrid (PLAXIS Output)
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
katai_add_test(test_plate_plastic katai::analysis)     # Plate M-N hinge (Mp/Np diamond, MMM section 18.3; section 10)
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
katai_add_test(test_axisym_gui katai::jobs katai::mesh katai::analysis katai::model)  # Axisymmetric GUI path: K0 + r-weighted gravity/traction + radial hoop
katai_add_test(test_gui_matrix katai::jobs katai::mesh katai::analysis katai::model)  # Systematic: structures x load x water x model x tri6/tri15 invariants
katai_add_test(test_safety_gui katai::jobs katai::mesh katai::analysis katai::model)  # Safety (phi-c reduction) GUI path: slope FoS + collapse mechanism
katai_add_test(test_k0_slope katai::jobs katai::mesh katai::analysis katai::model)    # K0 nil-step equilibrium on a sloping surface (horizontal = exact u0; a real redistribution on the slope)
katai_add_test(test_seepage_gui katai::jobs katai::mesh katai::analysis katai::model) # Groundwater flow GUI path: exact Darcy + Charny seepage face + Terzaghi coupling
katai_add_test(test_consolidation_gui katai::jobs katai::mesh katai::analysis katai::model) # Consolidation GUI path (PhaseType::Consolidation): 1D Terzaghi U-Tv + final settlement
katai_add_test(test_dynamic_gui katai::jobs katai::mesh katai::analysis katai::model) # Dynamic/seismic GUI path (PhaseType::Dynamic, D4b): 1D site-response resonance (4/pi A/w1^2 2xi) + frequency selectivity
katai_add_test(test_water_gui katai::jobs katai::mesh katai::analysis katai::model) # Water-history GUI path: FullyCoupled = Terzaghi + TransientFlow flow (W4 wiring)
katai_add_test(test_undrained_gui katai::jobs katai::mesh katai::analysis katai::model) # Undrained (A/B) GUI path: confined column -qH/M_u + effective stress (Skempton 1D)
katai_add_test(test_robustness katai::jobs katai::mesh katai::analysis katai::model) # Robustness: the calculation path never crashes or leaks NaN on hostile/degenerate input
katai_add_test(test_struct_gui katai::jobs katai::mesh katai::analysis katai::model)  # Structural M/Q/N diagrams from the GUI path (wiring: stations along the line, free-end M0, strut N)
katai_add_test(test_interface_gui katai::jobs katai::mesh katai::analysis katai::model)  # Standalone interface GUI path: embed in mesh + split + soil-soil Coulomb -> present != absent (tri6+tri15)
katai_add_test(test_mesh_refine katai::jobs katai::mesh katai::model) # Local mesh density (coarseness/sizing field): region + line refinement, quality, bit-for-bit legacy
katai_add_test(test_staged_gui katai::jobs katai::mesh katai::analysis katai::model)  # Multi-phase staged construction (SumMstage): fill settlement + excavation heave + phase-chain load (closed form)
katai_add_test(test_staged_struct_carry katai::jobs katai::mesh katai::analysis katai::model)  # Static chain structural state carry: NIL phase no-op IDENTITY (round-off; caught a 32% M shift on 2026-07-19 -> fixed)
katai_add_test(test_site_response_benchmark katai::jobs katai::mesh katai::analysis katai::model)  # BENCHMARK Wave-1: multi-layer site response vs a SHAKE-class layered viscoelastic transfer matrix (independent oracle, Track 7)
katai_add_test(test_compliant_base katai::jobs katai::mesh katai::analysis katai::model)  # Compliant (absorbing) base vs a layered solution with a radiation BC: uniform |T|~1, rigid-vs-compliant tooth, honest guards
katai_add_test(test_real_record katai::jobs katai::mesh katai::analysis katai::model)  # Real record: accelerogram input (bit-for-bit wiring) + El Centro 1940 NS (published PGA + spectral shape) + a product run
katai_add_test(test_plaxis_validation katai::mesh katai::analysis katai::geometry)  # Official PLAXIS Validation Manual cases: Giroud footing, Gibson E(z), Cox circular bearing, Davis-Booker c(z) (analytic | PLAXIS | KATAI)
katai_add_test(test_soft_soil katai::mesh katai::analysis katai::geometry)  # Soft Soil: stage 1 material point (ln law lambda*/kappa*, memory, K0NC pin, MC failure) + stage 2 FE wiring (integrate_point identity, frame independence, axisymmetry, oedometer BVP) (MMM section 10)
katai_add_test(test_soft_soil_gui katai::jobs katai::mesh katai::analysis katai::model)  # Soft Soil stage 3 GUI path: K0 + ss_initial_pp seeding, NC lambda* integral, OCR/POP memory, honest undrained refusal, dynamic p'/kappa* stiffness
katai_add_test(test_soft_soil_creep katai::mesh katai::analysis katai::geometry)  # Soft Soil Creep: stage 1 material point (mu*.ln(1+t/tau), OCR^-beta, closed-form relaxation, 1-day NC line, MC) + stage 2 FE (time into integrate_point, time-stepped relaxation BVP identity) (MMM section 11)
katai_add_test(test_input_audit katai::jobs katai::mesh katai::analysis katai::model)  # Input audit pins: OCR -> K0 (PLAXIS formula), HS+Undrained(B) refusal, plate Mp/Np refusal
if(TARGET test_real_record)
  target_compile_definitions(test_real_record PRIVATE KATAI_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
endif()
katai_add_test(test_design_gui katai::jobs katai::mesh katai::analysis katai::model)  # EC7 design approach through the GUI path (build_problem): DA3 ODF = FoS/gamma_M identity + None baseline (v0.3 B1.3)
katai_add_test(test_project_io katai::model katai::io)                # Project file JSON round trip (every field exact + honest failure modes)
katai_add_test(test_results_io katai::mesh katai::analysis katai::model katai::io)  # Binary result file: bitwise round trip + model-hash staleness + refusal of a corrupt file
katai_add_test(test_k2d_validator katai::model katai::io)             # .k2d input contract: every rule probed at its exact field path + reader forward-version notes
katai_add_test(test_input_corpus katai::jobs katai::mesh katai::analysis katai::model katai::io)  # Input corpus (roadmap 6.3): checked-in .k2d byte-identical to the programmatic build, validated, then SOLVED FROM THE FILE against its declared closed form
if(TARGET test_input_corpus)
  target_compile_definitions(test_input_corpus PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus"
    KATAI_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")   # El Centro record (KV-DYN-003)
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
katai_add_test(test_no_silent_drop katai::jobs katai::mesh katai::analysis katai::model katai::io)  # KV-DIA-001: no input may be discarded in silence -- one fixture per driver build loop, each perturbing ONE field of a corpus case, asserting the stable diagnostic code and that a clipped object equals the object drawn where the mesh took it
if(TARGET test_no_silent_drop)
  target_compile_definitions(test_no_silent_drop PRIVATE
    KATAI_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus")
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
add_executable(study_slope_convergence EXCLUDE_FROM_ALL study_slope_convergence.cpp)  # mesh convergence of the slope SRM -> docs/validation/slope-srm-convergence.md
target_link_libraries(study_slope_convergence PRIVATE katai::analysis katai::linsolve)
add_executable(study_design_ec7 EXCLUDE_FROM_ALL study_design_ec7.cpp)  # EC7 DA3 material factoring (ODF=FoS/gamma_M identity + Griffiths & Lane benchmark + PLAXIS algorithm identity)
target_link_libraries(study_design_ec7 PRIVATE katai::analysis katai::mesh katai::linsolve)
add_executable(study_wall_benchmark EXCLUDE_FROM_ALL study_wall_benchmark.cpp)  # PLAXIS sheet-pile clay benchmark
target_link_libraries(study_wall_benchmark PRIVATE katai::analysis katai::linsolve)
add_executable(study_footing_plaxis EXCLUDE_FROM_ALL study_footing_plaxis.cpp)  # PLAXIS Tutorial L1 circular footing
target_link_libraries(study_footing_plaxis PRIVATE katai::analysis katai::linsolve)
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

# Say out loud what this composition left out, so a shrunken test set is a stated
# fact rather than a surprise in the ctest total.
get_property(_katai_skipped GLOBAL PROPERTY KATAI_SKIPPED_TESTS)
if(_katai_skipped)
  list(LENGTH _katai_skipped _n)
  message(STATUS "katai tests: ${_n} skipped (their modules are not in this composition): ${_katai_skipped}")
endif()
