#pragma once
// WHAT A RUNNING PHASE IS DOING, for whoever is watching it. A front end gives Job a place to read
// from; solve_phases installs it on the solving thread for the length of the run, and the solver writes
// to it at every iteration (how far the iterate is from equilibrium, how many stress points are plastic)
// and at every step it accepts (the factor it reached and the largest movement, which is the phase's
// curve). A reader on another thread takes a copy.
//
// Nothing here reaches the answer: a note writes what the solver already has in hand. Without a scope
// (tests, the Python API, a plain solve call) a note is a load and a branch.

#include <array>
#include <mutex>
#include <vector>

namespace katai::analysis {

class LiveProgress {
public:
    struct Snapshot {
        int step = 0;            // steps accepted in this phase
        int iteration = 0;       // iterations of the step under way
        int plastic_points = 0;  // stress points plastic at the last iterate
        double residual = 0.0;   // force out of balance at the last iterate, over the phase's force scale
        double factor = 0.0;     // load factor, or strength factor in a Safety phase, of the last step
        double movement = 0.0;   // largest movement of a node in this phase, at the last step [m]
        std::vector<std::array<double, 2>> curve;   // {movement, factor} of every step accepted
    };

    void start_phase() {
        const std::lock_guard<std::mutex> lock(mu_);
        now_ = Snapshot{};
    }
    void iterate(int iteration, double residual, int plastic_points) {
        const std::lock_guard<std::mutex> lock(mu_);
        now_.iteration = iteration;
        now_.residual = residual;
        if (plastic_points > 0) now_.plastic_points = plastic_points;   // an iterate that kept its assembly counts none
    }
    void step(double factor, double movement) {
        const std::lock_guard<std::mutex> lock(mu_);
        ++now_.step;
        now_.factor = factor;
        now_.movement = movement;
        now_.curve.push_back({movement, factor});
    }
    Snapshot snapshot() const {
        const std::lock_guard<std::mutex> lock(mu_);
        return now_;
    }

private:
    mutable std::mutex mu_;
    Snapshot now_;
};

namespace detail {
inline LiveProgress*& live_progress() {
    thread_local LiveProgress* live = nullptr;
    return live;
}
}  // namespace detail

inline void note_iterate(int iteration, double residual, int plastic_points) {
    if (LiveProgress* live = detail::live_progress()) live->iterate(iteration, residual, plastic_points);
}
inline void note_step(double factor, double movement) {
    if (LiveProgress* live = detail::live_progress()) live->step(factor, movement);
}

// Installs `live` on the calling thread until the scope closes. Scopes nest.
class LiveScope {
public:
    explicit LiveScope(LiveProgress* live) : prev_(detail::live_progress()) { detail::live_progress() = live; }
    ~LiveScope() { detail::live_progress() = prev_; }
    LiveScope(const LiveScope&) = delete;
    LiveScope& operator=(const LiveScope&) = delete;

private:
    LiveProgress* prev_;
};

}  // namespace katai::analysis
