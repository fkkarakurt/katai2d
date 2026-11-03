#pragma once
// Cooperative cancellation INSIDE a phase. A front end asks a run to stop (Job::request_cancel);
// solve_phases installs that request on the solving thread for the length of the run, and the
// step, iteration and line-search loops of every solver poll it. A poll that finds the request
// set throws Cancelled. It unwinds the phase being solved, whose partial state is discarded and
// never reported, to solve_phases, which returns the phases that had completed.
//
// A poll between iterations is not enough: ONE assembly far from equilibrium can take half a
// minute (measured: the third excavation stage of a staged excavation in Hardening Soil sand, a
// line-search trial whose stress integration took 28 s), and a Cancel waited for it. The
// parallel assembly loops therefore poll per element too, through a CancelRegion: the workers do
// not carry the solving thread's request (it is thread_local), so the region takes it from the
// thread that opens it, and asks it one thread at a time, at most every 10 ms -- a request need
// not be safe to call from two threads at once -- and once it answers yes the other threads stop
// without asking again.
//
// Cancelled is deliberately NOT a std::exception: the handlers that turn a solver's exception
// into a failed phase ("Calculation failed -- ...") must not report a cancelled run as a
// failure. A catch-all on the path between a poll and solve_phases has to rethrow it.
//
// Without a scope (tests, the Python API, a plain solve call) a poll is a load and a branch and
// nothing can stop a solve.

#include <katai/materials/interrupt.hpp>

#include <atomic>
#include <chrono>
#include <functional>
#include <mutex>

namespace katai::analysis {

struct Cancelled {};

namespace detail {
inline const std::function<bool()>*& cancel_request() {
    thread_local const std::function<bool()>* request = nullptr;
    return request;
}
}  // namespace detail

inline void poll_cancel() {
    const auto* request = detail::cancel_request();
    if (request && (*request)()) throw Cancelled{};
}

// The request of the thread that opens a parallel region, polled from any thread of the region
// (see above). Without a scope on the opening thread, poll() is a load and a branch.
class CancelRegion {
public:
    CancelRegion() : request_(detail::cancel_request()) {}
    CancelRegion(const CancelRegion&) = delete;
    CancelRegion& operator=(const CancelRegion&) = delete;

    void poll() {
        if (!request_) return;
        if (stopped_.load(std::memory_order_relaxed)) throw Cancelled{};
        const long long now = std::chrono::duration_cast<std::chrono::microseconds>(
                                  std::chrono::steady_clock::now().time_since_epoch()).count();
        if (now < next_.load(std::memory_order_relaxed)) return;
        std::unique_lock<std::mutex> lock(mu_, std::try_to_lock);
        if (!lock.owns_lock() || stopped_.load(std::memory_order_relaxed)) return;
        next_.store(now + kIntervalUs, std::memory_order_relaxed);
        if ((*request_)()) {
            stopped_.store(true, std::memory_order_relaxed);
            throw Cancelled{};
        }
    }

    // Installs the region as the calling thread's interrupt hook (katai/materials/interrupt.hpp)
    // for the length of a piece of the region's work: a stress-point walk of a thousand substeps
    // polls between them. Seconds in one point were measured where the Hardening Soil walk ran to
    // its ceiling, and the poll per element could not see inside them.
    class Bind {
    public:
        explicit Bind(CancelRegion& region) : prev_(core::interrupt_hook()) {
            if (region.request_) core::interrupt_hook() = {&CancelRegion::hook, &region};
        }
        ~Bind() { core::interrupt_hook() = prev_; }
        Bind(const Bind&) = delete;
        Bind& operator=(const Bind&) = delete;

    private:
        core::InterruptHook prev_;
    };

private:
    static void hook(void* region) { static_cast<CancelRegion*>(region)->poll(); }

    static constexpr long long kIntervalUs = 10000;
    const std::function<bool()>* request_;
    std::atomic<bool> stopped_{false};
    std::atomic<long long> next_{0};
    std::mutex mu_;
};

// Installs `request` on the calling thread until the scope closes (an empty function installs
// nothing). Scopes nest; the inner one wins and the outer one comes back when it closes. The
// scope is also the thread's interrupt hook, so that the stress-point walks the solving thread
// runs outside a parallel region -- the cap calibration of a Hardening Soil before the first
// phase, seconds on a staged excavation with two of them -- stop with the run as well.
class CancelScope {
public:
    explicit CancelScope(const std::function<bool()>& request)
        : prev_(install(request)), bind_(region_) {}
    ~CancelScope() { detail::cancel_request() = prev_; }
    CancelScope(const CancelScope&) = delete;
    CancelScope& operator=(const CancelScope&) = delete;

private:
    static const std::function<bool()>* install(const std::function<bool()>& request) {
        const std::function<bool()>* prev = detail::cancel_request();
        if (request) detail::cancel_request() = &request;
        return prev;
    }

    const std::function<bool()>* prev_;
    CancelRegion region_;   // reads the request just installed
    CancelRegion::Bind bind_;
};

}  // namespace katai::analysis
