#pragma once
// The vocabulary the engine and the solver backends share for a failed linear solve.
//
// An analysis has to tell two situations apart, because one is recoverable and the
// other is not:
//
//   SingularSystem  The system has no usable solution. Either the factorization
//                   refused it, or an answer came back that does not satisfy it. At a
//                   limit load this is normal rather than exceptional: once a collapse
//                   mechanism forms, the tangent is rank-deficient along that
//                   mechanism and the out-of-balance force has a component along it,
//                   so no Newton direction exists. The analysis recovers by abandoning
//                   the increment and cutting it back, and the load level it last
//                   equilibrated is the answer being sought.
//
//   SolveError      The solve could not be attempted -- a matrix that is not square,
//                   an empty system, a solve before a factorization. There is nothing
//                   to recover from, and continuing would turn a fault into a reported
//                   collapse load: a plausible number produced by a bug, which is the
//                   worst thing this code could publish.
//
// Catching std::exception in the Newton loop would collapse that distinction, which is
// why the distinction is a type rather than a convention.
//
// Note which axis this splits on. It is *not* "did the factorization fail or did the
// verification fail" -- that is a property of the backend, not of the model. A
// singular tangent makes Eigen's factorization refuse outright, while PARDISO perturbs
// the pivots, completes, and produces an answer the residual check then rejects. Same
// model, same defect, two different routes to it, so both raise SingularSystem and the
// engine behaves identically on either backend.
//
// These live in katai::math because math sits below both katai::core and
// katai::linsolve, so neither has to depend on the other in order to name them.

#include <stdexcept>
#include <string>

namespace katai::math {

class SolveError : public std::runtime_error {
public:
    explicit SolveError(const std::string& what) : std::runtime_error(what) {}
};

// Deliberately derived from SolveError: a caller that does not care about the
// distinction, and every handler written before it existed, keeps working unchanged.
class SingularSystem : public SolveError {
public:
    explicit SingularSystem(const std::string& what) : SolveError(what) {}
};

}  // namespace katai::math
