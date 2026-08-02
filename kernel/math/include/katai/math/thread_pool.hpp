#pragma once
// Persistent thread pool + parallel_for — for the parallel phase of the FEM assembly.
//
// Why not OpenMP: MSVC OpenMP (vcomp) workers SPIN-WAIT after a parallel region ends and
// race for cores with the single-threaded PARDISO factorization that runs right after
// (measured: PARDISO ~2.4× slower). This pool's workers SLEEP on a condition_variable
// while waiting → zero cost to serial sections; there is also no external
// runtime/DLL dependency (the project's dependency-free philosophy).
//
// Determinism: parallel_for splits [0,n) into FIXED contiguous blocks per call; the
// caller works too. Do NOT accumulate shared floating point — each block must write only
// its own output (two-phase assembly: compute parallel, scatter sequential); then the
// result is bit-for-bit independent of the thread count.

#include <functional>

namespace katai::math {

// fn(begin, end): processes the contiguous subrange [begin,end). Blocks run concurrently
// on the pool workers + the calling thread; on return all have finished. If a block
// throws, the first exception is rethrown to the caller (the other blocks complete).
// For small n (or when called from inside a pool worker) it runs directly on the caller.
void parallel_for(int n, const std::function<void(int, int)>& fn);

// The pool's total parallel width (workers + caller). At least 1.
int parallel_threads();

} // namespace katai::math
