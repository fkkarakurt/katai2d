// CSR matrix + SparseMatrixBuilder unit test.
// CTest + assert style: 0 = passed, 1 = failed.
#include <katai/math/sparse_matrix.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::math::CsrMatrix;
using katai::math::Index;
using katai::math::Scalar;
using katai::math::SparseMatrixBuilder;

namespace {

int g_failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++g_failures;
    }
}

bool close(Scalar a, Scalar b) { return std::fabs(a - b) < 1e-12; }

// Summing of duplicate entries + in-row column sorting.
void test_build_sorts_and_sums_duplicates() {
    SparseMatrixBuilder builder(3, 3);
    builder.add_entry(0, 2, 1.0);
    builder.add_entry(0, 0, 2.0);
    builder.add_entry(0, 2, 3.0);  // (0,2) repeated → must sum to 1+3 = 4
    builder.add_entry(2, 1, 5.0);
    builder.add_entry(1, 1, 7.0);

    const CsrMatrix m = builder.build();

    check(m.rows == 3 && m.cols == 3, "dimensions");
    check(m.nonzeros() == 4, "nnz = 4 (one pair merged)");

    const std::vector<Index> expected_row_ptr = {0, 2, 3, 4};
    const std::vector<Index> expected_cols = {0, 2, 1, 1};
    const std::vector<Scalar> expected_vals = {2.0, 4.0, 7.0, 5.0};
    check(m.row_ptr == expected_row_ptr, "row_ptr");
    check(m.col_indices == expected_cols, "col_indices (sorted within rows)");
    check(m.values.size() == expected_vals.size(), "values size");
    for (std::size_t k = 0; k < expected_vals.size(); ++k)
        check(close(m.values[k], expected_vals[k]), "values content");
}

// y = A * x verification.
void test_multiply() {
    SparseMatrixBuilder builder(3, 3);
    builder.add_entry(0, 0, 2.0);
    builder.add_entry(0, 2, 4.0);
    builder.add_entry(1, 1, 7.0);
    builder.add_entry(2, 1, 5.0);
    const CsrMatrix m = builder.build();

    const std::vector<Scalar> x = {1.0, 1.0, 1.0};
    std::vector<Scalar> y(3, 0.0);
    m.multiply(x.data(), y.data());
    check(close(y[0], 6.0) && close(y[1], 7.0) && close(y[2], 5.0), "multiply");
}

// add_block: scattering a dense element block to global DOFs.
void test_add_block_scatter() {
    SparseMatrixBuilder builder(3, 3);
    const Index dof_map[2] = {1, 2};
    const Scalar block[4] = {10.0, 20.0,
                             30.0, 40.0};  // 2x2 row-major
    builder.add_block(dof_map, 2, block);
    const CsrMatrix m = builder.build();

    // Expected: (1,1)=10 (1,2)=20 (2,1)=30 (2,2)=40
    const std::vector<Index> expected_row_ptr = {0, 0, 2, 4};
    check(m.row_ptr == expected_row_ptr, "block row_ptr (row 0 empty)");
    check(m.nonzeros() == 4, "block nnz");
    std::vector<Scalar> y(3, 0.0);
    const std::vector<Scalar> x = {1.0, 1.0, 1.0};
    m.multiply(x.data(), y.data());
    check(close(y[1], 30.0) && close(y[2], 70.0), "block multiply");
}

// build_cached: the fast path (same pattern signature) must produce a matrix
// BIT-IDENTICAL to build() (summation order included); on a pattern change it must fall to
// a full build by itself.
void test_build_cached() {
    katai::math::CsrPatternCache cache;
    SparseMatrixBuilder b(3);
    auto fill = [&](Scalar s) {
        b.clear();
        // Duplicate entries + unsorted columns (mimics the FEM scatter).
        b.add_entry(0, 2, 3.0 * s); b.add_entry(0, 0, 1.0 * s); b.add_entry(0, 2, 0.5 * s);
        b.add_entry(1, 1, 2.0 * s); b.add_entry(2, 0, 4.0 * s); b.add_entry(1, 1, -1.0 * s);
        b.add_entry(2, 2, 5.0 * s);
    };
    fill(1.0);
    const CsrMatrix ref1 = b.build();
    const CsrMatrix& c1 = b.build_cached(cache);   // first: full build + map
    check(c1.row_ptr == ref1.row_ptr && c1.col_indices == ref1.col_indices &&
          c1.values == ref1.values, "build_cached first call identical to build");

    fill(7.5);                                      // same pattern, new values
    const CsrMatrix ref2 = b.build();
    const CsrMatrix& c2 = b.build_cached(cache);   // fast path (no sort)
    check(c2.row_ptr == ref2.row_ptr && c2.col_indices == ref2.col_indices &&
          c2.values == ref2.values, "build_cached fast path bit-for-bit");

    b.clear();                                      // PATTERN changed → fall to full build
    b.add_entry(0, 1, 9.0); b.add_entry(2, 2, 8.0);
    const CsrMatrix ref3 = b.build();
    const CsrMatrix& c3 = b.build_cached(cache);
    check(c3.row_ptr == ref3.row_ptr && c3.col_indices == ref3.col_indices &&
          c3.values == ref3.values, "build_cached full build on a pattern change");
}

} // namespace

int main() {
    test_build_sorts_and_sums_duplicates();
    test_multiply();
    test_add_block_scatter();
    test_build_cached();

    if (g_failures == 0) {
        std::printf("OK: sparse_matrix passed all checks\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
