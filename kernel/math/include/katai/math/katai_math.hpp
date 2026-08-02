#pragma once
// katai_math — geoteknik FEM math kernel (iskelet).
// Phase 0 / P0.1: version info + an Eigen smoke test only.

namespace katai::math {

// Module version string (proof of cross-module linking).
const char* version();

// A small SPD solve proving Eigen really compiles and runs.
// Solves a 2x2 system, returns the residual ||A x - b|| (≈ 0 expected).
double smoke_solve_residual();

} // namespace katai::math
