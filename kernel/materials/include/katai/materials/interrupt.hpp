#pragma once
// A hook a long stress-point integration calls between its substeps, so that a run can be stopped
// inside one point's walk rather than after it. The materials library knows nothing of runs: the
// analysis layer installs the hook on each thread of a parallel assembly (katai/analysis/cancel.hpp,
// CancelRegion) and the hook throws to stop the run. Where nothing is installed -- a material
// point driver, a test, a calibration -- a call is a load and a branch.

namespace katai::core {

struct InterruptHook {
    void (*fn)(void*) = nullptr;
    void* ctx = nullptr;
};

inline InterruptHook& interrupt_hook() {
    thread_local InterruptHook hook;
    return hook;
}

inline void poll_interrupt() {
    const InterruptHook& h = interrupt_hook();
    if (h.fn) h.fn(h.ctx);
}

}  // namespace katai::core
