#pragma once
// Sparse matrix types — for the FEM global stiffness matrix.
//
// Two parts:
//   * CsrMatrix          — zero-based CSR (Compressed Sparse Row), the final form
//                          handed to the solver.
//   * SparseMatrixBuilder — accumulates (row, column, value) triplets (COO), then
//                          compiles to CSR in one pass. FEM assembly contributes to
//                          the same (row,column) pair many times; those are summed
//                          during build().

#include <cstddef>
#include <vector>

#include <Eigen/Core>

namespace katai::math {

// The index type aligns with the MKL lp64 interface (32-bit). On a move to ilp64 this
// becomes std::int64_t and the static_assert in pardiso_solver.cpp fires.
using Index = int;
using Scalar = double;

// Zero-based CSR matrix. Column indices ascend within each row; duplicate entries were
// summed during compilation (deduplicated).
struct CsrMatrix {
    Index rows = 0;
    Index cols = 0;
    std::vector<Index> row_ptr;       // size: rows + 1 (row start offsets)
    std::vector<Index> col_indices;   // size: nnz
    std::vector<Scalar> values;       // size: nnz

    Index nonzeros() const { return static_cast<Index>(values.size()); }
    bool is_square() const { return rows == cols; }

    // y = A * x  (for verification/tests; x and y at least `rows`/`cols` long).
    void multiply(const Scalar* x, Scalar* y) const;
    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const;
};

// CSR pattern cache — for build_cached(). Newton iterations refill the global matrix
// with the SAME COO entry order every time (the assembly loops are deterministic); since
// the pattern is unchanged, sorting/deduplication happens once and later iterations only
// accumulate values through the slot map.
struct CsrPatternCache {
    CsrMatrix matrix;                 // the pattern + the last built values
    // Scan-order map: in the ORDER build() accumulates values (row → column-sorted),
    // which COO entry was added to which value slot. The fast path accumulates in the
    // same order → the floating-point summation order is IDENTICAL to a full build
    // (not even round-off deviates).
    std::vector<Index> scan_entry;    // scan position i → COO entry index k
    std::vector<Index> scan_slot;     // scan position i → matrix.values slot
    std::vector<Index> entry_rows;    // the cache's validity signature: are the entry
    std::vector<Index> entry_cols;    //  (row, column) arrays exactly the same?
};

// Helper that accumulates COO triplets and produces a CSR matrix.
// Filling is cheap (just push_back); the cost is in the one-time build().
class SparseMatrixBuilder {
public:
    SparseMatrixBuilder(Index rows, Index cols);
    explicit SparseMatrixBuilder(Index n) : SparseMatrixBuilder(n, n) {}

    // Reserve room for the expected total entry count (reduces reallocation).
    void reserve(std::size_t entry_count);

    // Adds a single contribution. Repeats of the same (row, col) are summed in build().
    void add_entry(Index row, Index col, Scalar value);

    // Scatters a dense element block through the global DOF mapping (FEM scatter).
    // dof_map: maps local index 0..count-1 to the global row/column.
    // block:   count x count, row-major.
    void add_block(const Index* dof_map, Index count, const Scalar* block_row_major);

    void clear();
    Index rows() const { return rows_; }
    Index cols() const { return cols_; }
    std::size_t entry_count() const { return values_.size(); }

    // Produces the CSR matrix by sorting the triplets and summing duplicates.
    CsrMatrix build() const;

    // Produces the SAME matrix as build(), but caches the pattern: if the cache is
    // populated and the current COO entries (as row,column arrays) match the previous
    // ones exactly, sorting/deduplication is skipped and only the values are accumulated
    // through the slot map (O(nnz), no sort). On a signature mismatch a full build runs
    // and the cache is refreshed — the result is identical to build() in every case. The
    // returned reference is cache.matrix.
    const CsrMatrix& build_cached(CsrPatternCache& cache) const;

private:
    Index rows_;
    Index cols_;
    std::vector<Index> rows_idx_;
    std::vector<Index> cols_idx_;
    std::vector<Scalar> values_;
};

} // namespace katai::math
