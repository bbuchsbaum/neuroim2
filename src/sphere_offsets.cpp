// Sphere-offset template and per-centre expansion.
//
// The set of voxel offsets making up a spherical neighbourhood depends only on
// (radius, spacing) -- not on the centre. local_sphere() in indexFuns.cpp
// recomputes it for every centre, which dominates exhaustive searchlights.
// sphere_offsets_cpp() builds the template once; sphere_at_cpp() then only
// translates and clips it.
//
// Equivalence with local_sphere():
//   * The membership test is written with the *identical* floating-point
//     expression, so the two agree bit-for-bit on boundary voxels.
//   * local_sphere() scans a cube of half-width ceil(radius/min(spacing)) on
//     every axis. Here each axis uses ceil(radius/spacing[axis]) instead. No
//     voxel is lost: an offset i with |i| > ceil(radius/spacing[axis]) has
//     |i*spacing[axis]| > radius and therefore fails the distance test anyway.
//     For anisotropic voxels this cuts the number of candidates by ~3x.
//   * Emission order is i (x) outer, j (y), k (z) inner -- matching
//     local_sphere(), which makes the rows lexicographically sorted by
//     (x, y, z) and lets spherical_roi() skip a re-sort.

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix sphere_offsets_cpp(double radius, NumericVector spacing) {
    if (spacing.size() < 3) {
        stop("sphere_offsets_cpp: 'spacing' must have at least 3 elements");
    }
    for (int a = 0; a < 3; ++a) {
        if (!(spacing[a] > 0) || !R_finite(spacing[a])) {
            stop("sphere_offsets_cpp: 'spacing' must be positive and finite");
        }
    }
    if (!R_finite(radius) || radius < 0) {
        stop("sphere_offsets_cpp: 'radius' must be non-negative and finite");
    }

    // Casting the ratio to int is undefined once it exceeds INT_MAX, and on
    // x86-64 it yields INT_MIN -- which is R's NA_integer_, so the template
    // would silently be full of NA offsets. Reject before the cast.
    double wd[3];
    for (int a = 0; a < 3; ++a) {
        wd[a] = std::ceil(radius / spacing[a]);
        if (!R_finite(wd[a]) || wd[a] > 2147483000.0) {
            stop("sphere_offsets_cpp: radius/spacing exceeds the representable "
                 "neighbourhood size on axis %d", a + 1);
        }
    }
    const int wx = (int) wd[0];
    const int wy = (int) wd[1];
    const int wz = (int) wd[2];

    std::vector<int> ox, oy, oz;
    // Typical searchlight neighbourhoods are a few hundred voxels; this is a
    // starting hint that avoids the first few reallocations, not a bound.
    ox.reserve(256); oy.reserve(256); oz.reserve(256);

    for (int i = -wx; i <= wx; i++) {
        for (int j = -wy; j <= wy; j++) {
            for (int k = -wz; k <= wz; k++) {
                // Same expression as local_sphere() -- do not "simplify" to
                // squared-distance form, that changes boundary rounding.
                double dist = sqrt(pow(i * spacing[0], 2) +
                                   pow(j * spacing[1], 2) +
                                   pow(k * spacing[2], 2));
                if (dist <= radius) {
                    ox.push_back(i); oy.push_back(j); oz.push_back(k);
                }
            }
        }
    }

    IntegerMatrix out(ox.size(), 3);
    for (size_t i = 0; i < ox.size(); i++) {
        out(i, 0) = ox[i];
        out(i, 1) = oy[i];
        out(i, 2) = oz[i];
    }
    return out;
}

// Translate a template to one centre and drop out-of-bounds voxels.
// `centre` and `dim` are 1-based; the result is 0-based when `base0` is TRUE,
// which is the convention local_sphere() used. Voxel coordinates are integers,
// so the result is an IntegerMatrix -- callers assembling a ROIVolWindow need
// integer coords and would otherwise pay for a copy to convert.
// [[Rcpp::export]]
IntegerMatrix sphere_at_cpp(IntegerMatrix off, IntegerVector centre, IntegerVector dim,
                            bool base0 = true) {
    if (centre.size() < 3 || dim.size() < 3) {
        stop("sphere_at_cpp: 'centre' and 'dim' must have at least 3 elements");
    }
    // Rcpp's operator() is unchecked, so a template with the wrong shape would
    // read past the end of the SEXP rather than error.
    if (off.ncol() != 3) {
        stop("sphere_at_cpp: 'off' must have exactly 3 columns");
    }
    const int cx = centre[0], cy = centre[1], cz = centre[2];
    const int d0 = dim[0], d1 = dim[1], d2 = dim[2];
    const int n = off.nrow();

    std::vector<int> keep;
    keep.reserve(n);
    for (int i = 0; i < n; i++) {
        const int x = cx + off(i, 0);
        if (x < 1 || x > d0) continue;
        const int y = cy + off(i, 1);
        if (y < 1 || y > d1) continue;
        const int z = cz + off(i, 2);
        if (z < 1 || z > d2) continue;
        keep.push_back(i);
    }

    const int shift = base0 ? 1 : 0;
    IntegerMatrix out(keep.size(), 3);
    for (size_t r = 0; r < keep.size(); r++) {
        const int i = keep[r];
        out(r, 0) = cx + off(i, 0) - shift;
        out(r, 1) = cy + off(i, 1) - shift;
        out(r, 2) = cz + off(i, 2) - shift;
    }
    return out;
}
