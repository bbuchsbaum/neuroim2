// Per-centre and batched expansion of a spherical neighbourhood into voxel
// coordinates.
//
// The searchlight iterators are lazy on purpose: a whole-brain radius-8 mm
// searchlight over a 290k-voxel mask is ~257 voxels per centre, which is about
// 0.9 GB of coordinates if materialised all at once. So the goal here is not to
// emit everything up front -- it is to make the *per-element* call cheap enough
// that the lazy iterator stops being dominated by R and S4 overhead.
//
// sphere_coords_cpp() does the whole per-centre job in one pass: translate the
// cached offset template, clip to the volume, and (optionally) drop voxels that
// are not in the mask. The callers hoist the template and the dimensions out of
// their closures, so the work per element is this call and nothing else.
//
// sphere_coords_batch_cpp() is the same loop over many centres, for the code
// paths that are eager by contract (spherical_roi_set(), searchlight(eager =
// TRUE)) and would otherwise call back into R once per centre.

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

using namespace Rcpp;

namespace {

// Shared inner loop.
//
// `vp` is the flattened volume, or NULL when no values are available. When
// `use_mask` is true a voxel is kept only if its value is present and non-zero
// -- NA and NaN are dropped, because "nonzero" cannot be established for them.
// That single rule is why filtering lives here rather than being reimplemented
// per caller: the R and C++ versions of it used to disagree about NA.
//
// `vv` collects the kept values (1.0 when no volume was supplied). Returns the
// 1-based row of the centre voxel, or NA_INTEGER when the centre is not kept.
inline int expand_one(const IntegerMatrix& off,
                      int cx, int cy, int cz,
                      int d0, int d1, int d2,
                      const double* vp, bool use_mask,
                      std::vector<int>& xs,
                      std::vector<int>& ys,
                      std::vector<int>& zs,
                      std::vector<double>* vv) {
    const int n = off.nrow();
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    int centre_row = NA_INTEGER;

    xs.clear(); ys.clear(); zs.clear();
    if (vv) vv->clear();

    for (int i = 0; i < n; i++) {
        const int ox = off(i, 0), oy = off(i, 1), oz = off(i, 2);
        const int x = cx + ox;
        if (x < 1 || x > d0) continue;
        const int y = cy + oy;
        if (y < 1 || y > d1) continue;
        const int z = cz + oz;
        if (z < 1 || z > d2) continue;

        double v = 1.0;
        if (vp) {
            const R_xlen_t lin = (R_xlen_t)(x - 1) + (R_xlen_t)(y - 1) * d0
                                 + (R_xlen_t)(z - 1) * slice;
            v = vp[lin];
            if (use_mask && (ISNAN(v) || v == 0.0)) continue;
        }

        if (ox == 0 && oy == 0 && oz == 0) {
            centre_row = (int) xs.size() + 1;
        }
        xs.push_back(x); ys.push_back(y); zs.push_back(z);
        if (vv) vv->push_back(v);
    }
    return centre_row;
}

inline IntegerMatrix to_matrix(const std::vector<int>& xs,
                               const std::vector<int>& ys,
                               const std::vector<int>& zs) {
    IntegerMatrix out(xs.size(), 3);
    for (size_t r = 0; r < xs.size(); r++) {
        out(r, 0) = xs[r];
        out(r, 1) = ys[r];
        out(r, 2) = zs[r];
    }
    return out;
}

inline void check_args(const IntegerMatrix& off, const IntegerVector& dim) {
    if (off.ncol() != 3) {
        stop("searchlight core: 'off' must have exactly 3 columns");
    }
    if (dim.size() < 3) {
        stop("searchlight core: 'dim' must have at least 3 elements");
    }
    // NA_INTEGER is INT_MIN, so it is caught here too.
    if (dim[0] <= 0 || dim[1] <= 0 || dim[2] <= 0) {
        stop("searchlight core: 'dim' must be positive");
    }
}

// Number of voxels implied by `dim`, computed with division-based checks so the
// product can never wrap. It matters: dim() on a NeuroVol reports the NeuroSpace
// slot, which R does not keep in sync with the data array, so a bogus dim can
// reach here. Forming d0*d1*d2 first would overflow R_xlen_t, wrap negative,
// slip past a `length < nvox` test and index far outside the buffer.
inline R_xlen_t checked_nvox(int d0, int d1, int d2) {
    const R_xlen_t lim = R_XLEN_T_MAX;
    if ((R_xlen_t) d1 > lim / (R_xlen_t) d0) {
        stop("searchlight core: 'dim' implies more voxels than are representable");
    }
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    if ((R_xlen_t) d2 > lim / slice) {
        stop("searchlight core: 'dim' implies more voxels than are representable");
    }
    return slice * (R_xlen_t) d2;
}

// The value buffer must either be absent (length 0, meaning "no volume
// supplied") or cover the whole volume. A short buffer is a programming error,
// not something to silently skip over.
inline void check_buffer(R_xlen_t len, R_xlen_t nvox, const char* what) {
    if (len > 0 && len < nvox) {
        stop("searchlight core: '%s' is shorter than prod(dim)", what);
    }
}

} // anonymous namespace

// One centre. Returns an n x 3 integer matrix of 1-based voxel coordinates,
// ordered lexicographically by (x, y, z) -- the template's own order, which
// clipping and filtering preserve.
// [[Rcpp::export]]
IntegerMatrix sphere_coords_cpp(IntegerMatrix off, IntegerVector centre,
                                IntegerVector dim, NumericVector vals,
                                bool use_mask) {
    check_args(off, dim);
    if (centre.size() < 3) {
        stop("searchlight core: 'centre' must have at least 3 elements");
    }
    check_buffer(vals.size(), checked_nvox(dim[0], dim[1], dim[2]), "vals");
    const double* vp = vals.size() > 0 ? REAL(vals) : NULL;

    std::vector<int> xs, ys, zs;
    xs.reserve(off.nrow()); ys.reserve(off.nrow()); zs.reserve(off.nrow());
    expand_one(off, centre[0], centre[1], centre[2], dim[0], dim[1], dim[2],
               vp, use_mask, xs, ys, zs, NULL);
    return to_matrix(xs, ys, zs);
}

// Many centres in one pass. `centres` is an m x 3 integer matrix of 1-based
// coordinates. Returns a list of coordinate matrices plus the 1-based row of
// each centre within its own matrix (NA where the centre did not survive
// filtering).
// [[Rcpp::export]]
List sphere_coords_batch_cpp(IntegerMatrix off, IntegerMatrix centres,
                             IntegerVector dim, NumericVector vals,
                             bool use_mask) {
    check_args(off, dim);
    if (centres.ncol() != 3) {
        stop("searchlight core: 'centres' must have exactly 3 columns");
    }

    const int m = centres.nrow();
    const int d0 = dim[0], d1 = dim[1], d2 = dim[2];
    check_buffer(vals.size(), checked_nvox(d0, d1, d2), "vals");
    const double* vp = vals.size() > 0 ? REAL(vals) : NULL;

    List coords(m);
    IntegerVector centre_row(m);

    std::vector<int> xs, ys, zs;
    xs.reserve(off.nrow()); ys.reserve(off.nrow()); zs.reserve(off.nrow());

    for (int c = 0; c < m; c++) {
        centre_row[c] = expand_one(off, centres(c, 0), centres(c, 1), centres(c, 2),
                                   d0, d1, d2, vp, use_mask, xs, ys, zs, NULL);
        coords[c] = to_matrix(xs, ys, zs);
    }

    return List::create(_["coords"] = coords, _["center_row"] = centre_row);
}

// One centre, everything a ROIVolWindow needs. Folds the value gather and the
// centre-row search into the same pass as the coordinate expansion, so the R
// side does nothing but hand the pieces to the constructor.
//
// `vals` is the flattened volume; length 0 means "emit indicator 1s".
// [[Rcpp::export]]
List sphere_roi_at_cpp(IntegerMatrix off, IntegerVector centre, IntegerVector dim,
                       NumericVector vals, bool use_mask) {
    check_args(off, dim);
    if (centre.size() < 3) {
        stop("searchlight core: 'centre' must have at least 3 elements");
    }

    const int cx = centre[0], cy = centre[1], cz = centre[2];
    const int d0 = dim[0], d1 = dim[1], d2 = dim[2];
    const R_xlen_t nvox = checked_nvox(d0, d1, d2);
    check_buffer(vals.size(), nvox, "vals");
    const double* vp = vals.size() > 0 ? REAL(vals) : NULL;

    std::vector<int> xs, ys, zs;
    std::vector<double> vv;
    xs.reserve(off.nrow()); ys.reserve(off.nrow()); zs.reserve(off.nrow());
    vv.reserve(off.nrow());

    const int centre_row = expand_one(off, cx, cy, cz, d0, d1, d2,
                                      vp, use_mask, xs, ys, zs, &vv);

    // Parent index can exceed INT_MAX for volumes above 2^31 voxels; report NA
    // rather than narrowing, which is what the R builders do.
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    const R_xlen_t parent = (R_xlen_t)(cz - 1) * slice + (R_xlen_t)(cy - 1) * d0 + cx;
    const int parent_i = (parent > INT_MAX) ? NA_INTEGER : (int) parent;

    return List::create(
        _["coords"]       = to_matrix(xs, ys, zs),
        _["values"]       = NumericVector(vv.begin(), vv.end()),
        _["center_row"]   = centre_row,
        _["parent_index"] = parent_i);
}
