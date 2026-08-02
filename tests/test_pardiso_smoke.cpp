// PARDISO link + correctness smoke test (Decision D6).
//
// Goal: prove not just that oneMKL compiles and links, but that it solves a real
// sparse system CORRECTLY. 5x5 SPD tridiagonal matrix:
//     A = tridiag(-1, 2, -1)
// b was chosen so the solution is x = [1,1,1,1,1] (b = A·1 = [1,0,0,0,1]).
// PARDISO symmetric (mtype=2) expects upper-triangle CSR, one-based indices.

#include <mkl.h>

#include <cmath>
#include <cstdio>

int main() {
    MKLVersion ver;
    mkl_get_version(&ver);
    std::printf("oneMKL %d.%d.%d\n", ver.MajorVersion, ver.MinorVersion,
                ver.UpdateVersion);

    constexpr MKL_INT n = 5;
    // CSR (upper triangle only, symmetric), one-based:
    MKL_INT ia[n + 1] = {1, 3, 5, 7, 9, 10};
    MKL_INT ja[9]     = {1, 2, 2, 3, 3, 4, 4, 5, 5};
    double  a[9]      = {2, -1, 2, -1, 2, -1, 2, -1, 2};
    double  b[n]      = {1, 0, 0, 0, 1};
    double  x[n]      = {0};

    void*   pt[64]    = {0};
    MKL_INT iparm[64] = {0};
    MKL_INT mtype     = 2;            // real symmetric positive definite
    pardisoinit(pt, &mtype, iparm);
    iparm[34] = 0;                    // one-based (Fortran) indexing

    MKL_INT maxfct = 1, mnum = 1, nrhs = 1, msglvl = 0, error = 0;
    MKL_INT idum = 0;
    double  ddum = 0.0;

    // Phase 13: analysis + numerical factorization + solve (one call).
    MKL_INT phase = 13;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, b, x, &error);
    if (error != 0) {
        std::printf("PARDISO solve error: %d\n", static_cast<int>(error));
        return 2;
    }

    // Faz -1: dahili bellek serbest.
    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    double maxerr = 0.0;
    for (int i = 0; i < n; ++i)
        maxerr = std::fmax(maxerr, std::fabs(x[i] - 1.0));
    std::printf("PARDISO maxerr = %.3e\n", maxerr);

    return (maxerr < 1e-9) ? 0 : 1;
}
