// Single-pass gathers for series().
//
// Both replace R-level loops that allocate a fresh index vector (or a whole
// zero matrix plus a scatter) on every call. Callers must only take these
// paths when the underlying storage is a plain double vector/matrix -- the
// R-side methods check that and fall back otherwise.
//
// NOTE ON THE SEXP ARGUMENTS
// The data arguments are taken as bare SEXP, not Rcpp::NumericVector, and
// deliberately so. `x@.Data` on an S4 object that extends "array" yields a
// REALSXP that Rcpp's Vector<REALSXP> constructor duplicates -- for a
// 48x56x40x200 volume that is a 164 MB memcpy on *every* series() call, which
// would make the "fast" path far slower than the R loop it replaces. Taking
// SEXP and reading REAL() after an explicit type check is guaranteed
// copy-free.

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

using namespace Rcpp;

// Gather time-series from a dense [x, y, z, t] volume.
//
// `dim` is the 4-D dimension vector, `coords` an n x 3 matrix of 1-based voxel
// coordinates (already bounds-checked on the R side). Returns [t x n].
// [[Rcpp::export]]
NumericMatrix series_gather_dense(SEXP data, IntegerVector dim, IntegerMatrix coords) {
    if (TYPEOF(data) != REALSXP) {
        stop("series_gather_dense: 'data' must be a double vector");
    }
    if (dim.size() < 4) stop("series_gather_dense: 'dim' must have 4 elements");
    if (coords.ncol() != 3) stop("series_gather_dense: 'coords' must have 3 columns");

    const R_xlen_t d0 = dim[0], d1 = dim[1], d2 = dim[2];
    const R_xlen_t nt = dim[3];
    const int n = coords.nrow();

    // NA_INTEGER is INT_MIN, so it is caught here too.
    if (d0 <= 0 || d1 <= 0 || d2 <= 0 || nt < 0) {
        stop("series_gather_dense: 'dim' must be positive");
    }

    // dim() on a NeuroVec reports the NeuroSpace slot, which R does not keep in
    // sync with the length of the underlying array -- so `dim` can claim far
    // more elements than `data` holds. Build the strides with division-based
    // checks: forming d0*d1*d2*nt first can overflow R_xlen_t and wrap negative,
    // which would let a bogus dim slip past the guard and read out of bounds.
    const R_xlen_t xlen = Rf_xlength(data);
    if (d1 > xlen / d0) {
        stop("series_gather_dense: 'data' is shorter than prod(dim)");
    }
    const R_xlen_t slice = d0 * d1;
    if (d2 > xlen / slice) {
        stop("series_gather_dense: 'data' is shorter than prod(dim)");
    }
    const R_xlen_t nels = slice * d2;
    if (nt > 0 && nels > xlen / nt) {
        stop("series_gather_dense: 'data' is shorter than prod(dim)");
    }

    NumericMatrix out(nt, n);
    const double* dp = REAL(data);

    for (int c = 0; c < n; c++) {
        const R_xlen_t x = coords(c, 0), y = coords(c, 1), z = coords(c, 2);
        // Defensive: the R side validates, but never index out of bounds.
        if (x < 1 || x > d0 || y < 1 || y > d1 || z < 1 || z > d2) {
            stop("series_gather_dense: coordinate out of bounds");
        }
        const R_xlen_t lin = (x - 1) + (y - 1) * d0 + (z - 1) * slice;
        double* op = &out(0, c);
        for (R_xlen_t t = 0; t < nt; t++) {
            op[t] = dp[lin + t * nels];
        }
    }
    return out;
}

// Gather columns from a sparse [time x voxel] matrix.
//
// `mapped` holds, for each requested voxel, its column in `data` (1-based) or 0
// when the voxel lies outside the mask. Columns for 0 entries stay zero, which
// is what the R implementation produced by pre-filling a zero matrix.
// Returns [time x length(mapped)].
// [[Rcpp::export]]
NumericMatrix series_gather_sparse(SEXP data, IntegerVector mapped) {
    if (TYPEOF(data) != REALSXP) {
        stop("series_gather_sparse: 'data' must be a double matrix");
    }
    SEXP dimAttr = Rf_getAttrib(data, R_DimSymbol);
    if (Rf_isNull(dimAttr) || Rf_length(dimAttr) != 2) {
        stop("series_gather_sparse: 'data' must be a matrix");
    }
    const R_xlen_t nt = INTEGER(dimAttr)[0];
    const R_xlen_t nk = INTEGER(dimAttr)[1];
    const int n = mapped.size();

    if (nt < 0 || nk < 0 || nt * nk > Rf_xlength(data)) {
        stop("series_gather_sparse: 'data' dim does not match its length");
    }

    NumericMatrix out(nt, n);
    const double* dp = REAL(data);

    for (int c = 0; c < n; c++) {
        const int m = mapped[c];
        if (m == NA_INTEGER || m <= 0) continue;   // leave the column zero
        if (m > nk) stop("series_gather_sparse: mapped index exceeds ncol(data)");
        const double* src = dp + (R_xlen_t)(m - 1) * nt;
        double* dst = &out(0, c);
        for (R_xlen_t t = 0; t < nt; t++) dst[t] = src[t];
    }
    return out;
}
