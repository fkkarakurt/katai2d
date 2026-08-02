#include <katai/math/sparse_matrix.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>

namespace katai::math {

// --- CsrMatrix -------------------------------------------------------------

void CsrMatrix::multiply(const Scalar* x, Scalar* y) const {
    for (Index r = 0; r < rows; ++r) {
        Scalar sum = Scalar{0};
        for (Index k = row_ptr[r]; k < row_ptr[r + 1]; ++k)
            sum += values[k] * x[col_indices[k]];
        y[r] = sum;
    }
}

Eigen::VectorXd CsrMatrix::operator*(const Eigen::VectorXd& x) const {
    assert(x.size() == cols);
    Eigen::VectorXd y(rows);
    multiply(x.data(), y.data());
    return y;
}

// --- SparseMatrixBuilder ---------------------------------------------------

SparseMatrixBuilder::SparseMatrixBuilder(Index rows, Index cols)
    : rows_(rows), cols_(cols) {
    assert(rows >= 0 && cols >= 0);
}

void SparseMatrixBuilder::reserve(std::size_t entry_count) {
    rows_idx_.reserve(entry_count);
    cols_idx_.reserve(entry_count);
    values_.reserve(entry_count);
}

void SparseMatrixBuilder::add_entry(Index row, Index col, Scalar value) {
    assert(row >= 0 && row < rows_ && col >= 0 && col < cols_);
    rows_idx_.push_back(row);
    cols_idx_.push_back(col);
    values_.push_back(value);
}

void SparseMatrixBuilder::add_block(const Index* dof_map, Index count,
                                    const Scalar* block_row_major) {
    for (Index i = 0; i < count; ++i)
        for (Index j = 0; j < count; ++j)
            add_entry(dof_map[i], dof_map[j], block_row_major[i * count + j]);
}

void SparseMatrixBuilder::clear() {
    rows_idx_.clear();
    cols_idx_.clear();
    values_.clear();
}

namespace {

// The shared body of build(). If cache is given, also produces the scan-order map:
// in the order build accumulates values (row → column-sorted), the COO entry index +
// the target value slot. When the fast path accumulates in this order, the
// floating-point summation order stays identical to a full build.
CsrMatrix build_csr(Index rows, Index cols,
                    const std::vector<Index>& rows_idx,
                    const std::vector<Index>& cols_idx,
                    const std::vector<Scalar>& values,
                    CsrPatternCache* cache) {
    const std::size_t nnz_in = values.size();

    CsrMatrix out;
    out.rows = rows;
    out.cols = cols;
    out.row_ptr.assign(static_cast<std::size_t>(rows) + 1, 0);

    // 1) Count entries per row → row_ptr (prefix sum).
    for (std::size_t k = 0; k < nnz_in; ++k)
        ++out.row_ptr[rows_idx[k] + 1];
    for (Index r = 0; r < rows; ++r)
        out.row_ptr[r + 1] += out.row_ptr[r];

    // 2) Scatter the entries into row buckets (not yet column-sorted).
    std::vector<Index> scratch_col(nnz_in);
    std::vector<Scalar> scratch_val(nnz_in);
    std::vector<Index> bucket_entry;            // bucket position → COO entry index
    if (cache) bucket_entry.resize(nnz_in);
    std::vector<Index> next = out.row_ptr;  // write cursor per row
    for (std::size_t k = 0; k < nnz_in; ++k) {
        const Index r = rows_idx[k];
        const Index dst = next[r]++;
        scratch_col[dst] = cols_idx[k];
        scratch_val[dst] = values[k];
        if (cache) bucket_entry[dst] = static_cast<Index>(k);
    }

    // 3) Sort each row by column and compress, summing duplicates.
    out.col_indices.reserve(nnz_in);
    out.values.reserve(nnz_in);
    std::vector<Index> order;  // in-row sorting permutation
    std::vector<Index> compact_row_ptr(static_cast<std::size_t>(rows) + 1, 0);
    if (cache) {
        cache->scan_entry.clear(); cache->scan_entry.reserve(nnz_in);
        cache->scan_slot.clear();  cache->scan_slot.reserve(nnz_in);
    }

    for (Index r = 0; r < rows; ++r) {
        const Index begin = out.row_ptr[r];
        const Index end = out.row_ptr[r + 1];
        const Index n = end - begin;

        order.resize(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](Index a, Index b) {
            return scratch_col[begin + a] < scratch_col[begin + b];
        });

        // Walk in ascending column order; sum equal columns.
        Index last_col = -1;
        for (Index i = 0; i < n; ++i) {
            const Index src = begin + order[i];
            const Index col = scratch_col[src];
            const Scalar val = scratch_val[src];
            if (col == last_col) {
                out.values.back() += val;
            } else {
                out.col_indices.push_back(col);
                out.values.push_back(val);
                last_col = col;
            }
            if (cache) {
                cache->scan_entry.push_back(bucket_entry[src]);
                cache->scan_slot.push_back(static_cast<Index>(out.values.size()) - 1);
            }
        }
        compact_row_ptr[r + 1] = static_cast<Index>(out.values.size());
    }

    out.row_ptr = std::move(compact_row_ptr);
    return out;
}

} // namespace

CsrMatrix SparseMatrixBuilder::build() const {
    return build_csr(rows_, cols_, rows_idx_, cols_idx_, values_, nullptr);
}

const CsrMatrix& SparseMatrixBuilder::build_cached(CsrPatternCache& cache) const {
    // Signature check: are the entry (row, column) arrays EXACTLY the previous ones?
    // An O(nnz) integer comparison — far cheaper than sorting, and it makes the cache
    // unconditionally correct (does not lean on the deterministic-assembly assumption).
    const bool same = cache.scan_entry.size() == values_.size() &&
                      cache.entry_rows == rows_idx_ && cache.entry_cols == cols_idx_;
    if (same) {
        std::fill(cache.matrix.values.begin(), cache.matrix.values.end(), Scalar{0});
        // Accumulate in scan order = exactly build()'s summation order.
        for (std::size_t i = 0; i < cache.scan_entry.size(); ++i)
            cache.matrix.values[cache.scan_slot[i]] += values_[cache.scan_entry[i]];
        return cache.matrix;
    }
    cache.matrix = build_csr(rows_, cols_, rows_idx_, cols_idx_, values_, &cache);
    cache.entry_rows = rows_idx_;
    cache.entry_cols = cols_idx_;
    return cache.matrix;
}

} // namespace katai::math
