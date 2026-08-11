// Separable Gaussian blur.
//
// gaussian_blur_cpp() in indexFuns.cpp applies the full (2w+1)^3 kernel at every
// mask voxel. A Gaussian is a tensor product, so the same result is available in
// 3(2w+1) taps via three 1-D passes -- 15x fewer multiplies at window = 3.
//
// The masked default (normalize = TRUE) is separable too, as a ratio of two
// separable convolutions. The dense code accumulates
//
//     acc  = sum over in-bounds, in-mask u of  w(u-v) * x[u]
//     wsum = sum over in-bounds, in-mask u of  w(u-v)
//
// and returns acc/wsum. Writing m for the mask indicator and y for x masked to
// zero outside it, acc is (w * y)[v] and wsum is (w * m)[v], both zero-padded
// convolutions of a separable kernel. So out = blur(y) / blur(m).
//
// y is built with an explicit branch rather than as x*m: out-of-mask voxels are
// documented to be allowed to hold NaN, and 0 * NaN is NaN, which would leak
// missingness into every in-mask neighbour.
//
// Kernel normalisation matches gaussian_weights(): the 3-D kernel there is a
// product of three 1-D factors divided by the total, and because the product is
// separable that total factorises, so normalising each 1-D factor by its own sum
// gives the identical normalised 3-D kernel.
//
// Cost: three passes over the *whole volume*, versus the dense kernel's one pass
// over the mask only. Whether that is a win depends on both the window and the
// mask fraction, so the caller decides -- see .gaussian_blur_engine() in
// R/spat_filter.R. Working memory is 2 vectors of prod(dim) doubles, or 3 plus a
// byte-per-voxel mask when normalising, which the dense path does not need.

// [[Rcpp::depends(RcppParallel)]]
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#include <RcppParallel.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif
#include <memory>
#include <climits>
#include <cstring>
#include <vector>

using namespace Rcpp;

namespace {

// One axis of gaussian_weights()'s 1-D factor, normalised to sum 1. The
// sqrt(2*pi*sigma) constant in that function is per-axis and cancels here.
std::vector<double> axis_kernel(int window, double sigma, double spacing) {
    if (window > (INT_MAX - 1) / 2) stop("gaussian_blur_sep_cpp: 'window' is too large");
    const int sz = 2 * window + 1;
    std::vector<double> k(sz);
    const double denom = 2.0 * sigma * sigma;
    double tot = 0.0;
    for (int i = -window; i <= window; i++) {
        const double v = std::exp(-std::pow(i * spacing, 2) / denom);
        k[i + window] = v;
        tot += v;
    }
    for (int i = 0; i < sz; i++) k[i] /= tot;
    return k;
}

// First pass along x, reading the masked source without materialising it.
//
//   in == nullptr  ->  source is the mask indicator m
//   in != nullptr  ->  source is y, i.e. in[u] where m[u], else 0
//
// The branch is what keeps out-of-mask NaN from leaking: y is *not* in*m,
// because 0 * NaN is NaN.
void pass_x_masked(const double* in, const char* m, double* out,
                   int d0, int d1, int d2,
                   const std::vector<double>& k, int w) {
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    for (int z = 0; z < d2; z++) {
        for (int y = 0; y < d1; y++) {
            const R_xlen_t base = (R_xlen_t) y * d0 + (R_xlen_t) z * slice;
            for (int x = 0; x < d0; x++) {
                double acc = 0.0;
                const int lo = std::max(-w, -x), hi = std::min(w, d0 - 1 - x);
                for (int i = lo; i <= hi; i++) {
                    const R_xlen_t u = base + x + i;
                    if (!m[u]) continue;
                    acc += k[i + w] * (in ? in[u] : 1.0);
                }
                out[base + x] = acc;
            }
        }
    }
}

// in -> out along x. Zero padding outside the volume.
void pass_x(const double* in, double* out, int d0, int d1, int d2,
            const std::vector<double>& k, int w) {
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    for (int z = 0; z < d2; z++) {
        for (int y = 0; y < d1; y++) {
            const R_xlen_t base = (R_xlen_t) y * d0 + (R_xlen_t) z * slice;
            for (int x = 0; x < d0; x++) {
                double acc = 0.0;
                const int lo = std::max(-w, -x), hi = std::min(w, d0 - 1 - x);
                for (int i = lo; i <= hi; i++) acc += k[i + w] * in[base + x + i];
                out[base + x] = acc;
            }
        }
    }
}

// The y and z passes accumulate over whole contiguous runs rather than one
// voxel at a time. Walking a column with stride d0 (or d0*d1) touches a fresh
// cache line on every iteration; adding tap i into an entire x-row at once
// reads and writes contiguous memory and vectorises. The additions for a given
// output element still happen in the same order (i ascending), so the result is
// bit-identical to the per-voxel form this replaces.
void pass_y(const double* in, double* out, int d0, int d1, int d2,
            const std::vector<double>& k, int w) {
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    for (int z = 0; z < d2; z++) {
        const double* inz = in + (R_xlen_t) z * slice;
        double* outz = out + (R_xlen_t) z * slice;
        for (int y = 0; y < d1; y++) {
            double* orow = outz + (R_xlen_t) y * d0;
            const int lo = std::max(-w, -y), hi = std::min(w, d1 - 1 - y);
            const double* irow0 = inz + (R_xlen_t)(y + lo) * d0;
            const double k0 = k[lo + w];
            for (int x = 0; x < d0; x++) orow[x] = k0 * irow0[x];
            for (int i = lo + 1; i <= hi; i++) {
                const double kw = k[i + w];
                const double* irow = inz + (R_xlen_t)(y + i) * d0;
                for (int x = 0; x < d0; x++) orow[x] += kw * irow[x];
            }
        }
    }
}

void pass_z(const double* in, double* out, int d0, int d1, int d2,
            const std::vector<double>& k, int w) {
    const R_xlen_t slice = (R_xlen_t) d0 * d1;
    for (int z = 0; z < d2; z++) {
        double* oslice = out + (R_xlen_t) z * slice;
        const int lo = std::max(-w, -z), hi = std::min(w, d2 - 1 - z);
        const double* islice0 = in + (R_xlen_t)(z + lo) * slice;
        const double k0 = k[lo + w];
        for (R_xlen_t p = 0; p < slice; p++) oslice[p] = k0 * islice0[p];
        for (int i = lo + 1; i <= hi; i++) {
            const double kw = k[i + w];
            const double* islice = in + (R_xlen_t)(z + i) * slice;
            for (R_xlen_t p = 0; p < slice; p++) oslice[p] += kw * islice[p];
        }
    }
}

// blur(m) along one axis when every voxel is in the mask: the zero-padded
// convolution of the constant 1 is separable in closed form, and the factor at
// coordinate i is just the kernel weight that stays in bounds there.
std::vector<double> edge_factor(int d, const std::vector<double>& k, int window) {
    std::vector<double> c(d);
    for (int i = 0; i < d; i++) {
        double s = 0.0;
        const int lo = std::max(-window, -i), hi = std::min(window, d - 1 - i);
        for (int t = lo; t <= hi; t++) s += k[t + window];
        c[i] = s;
    }
    return c;
}

// Three passes, ping-ponging between `a` and `b`. Reads `src` on the first pass
// so the caller need not copy it. Result lands in `a` (3 passes: src->b, b->a,
// a->b ... so it lands in b; see the ordering below).
// Returns a pointer to whichever buffer holds the result.
double* separable3(const double* src, double* a, double* b,
                   int d0, int d1, int d2,
                   const std::vector<double>& kx,
                   const std::vector<double>& ky,
                   const std::vector<double>& kz, int w) {
    pass_x(src, a, d0, d1, d2, kx, w);
    pass_y(a,   b, d0, d1, d2, ky, w);
    pass_z(b,   a, d0, d1, d2, kz, w);
    return a;
}

} // anonymous namespace

// Separable equivalent of gaussian_blur_cpp(). Same arguments, same return: a
// numeric vector carrying `dim`, blurred at `mask_idx` and zero elsewhere.
// [[Rcpp::export]]
NumericVector gaussian_blur_sep_cpp(NumericVector arr, IntegerVector mask_idx,
                                    int window, double sigma,
                                    NumericVector spacing, bool normalize = true,
                                    bool full_mask = false) {
    // Read the attribute as a bare SEXP: converting a missing dim straight to
    // IntegerVector throws Rcpp's "Not compatible with requested type" instead
    // of saying what is actually wrong.
    SEXP dim_attr = arr.attr("dim");
    if (Rf_isNull(dim_attr) || Rf_xlength(dim_attr) != 3)
        stop("gaussian_blur_sep_cpp: 'arr' must have a 3-D dim attribute");
    IntegerVector dims = dim_attr;
    const int d0 = dims[0], d1 = dims[1], d2 = dims[2];
    if (d0 <= 0 || d1 <= 0 || d2 <= 0) stop("gaussian_blur_sep_cpp: 'dim' must be positive");
    if (window < 0) stop("gaussian_blur_sep_cpp: 'window' must be non-negative");
    if (!(sigma > 0) || !R_finite(sigma)) stop("gaussian_blur_sep_cpp: 'sigma' must be positive");
    if (spacing.size() < 3) stop("gaussian_blur_sep_cpp: 'spacing' must have 3 elements");

    const R_xlen_t nvox = (R_xlen_t) d0 * d1 * d2;
    if (nvox > arr.size()) stop("gaussian_blur_sep_cpp: 'arr' is shorter than prod(dim)");

    const std::vector<double> kx = axis_kernel(window, sigma, spacing[0]);
    const std::vector<double> ky = axis_kernel(window, sigma, spacing[1]);
    const std::vector<double> kz = axis_kernel(window, sigma, spacing[2]);

    const R_xlen_t nm = mask_idx.size();
    const double* ap = REAL(arr);

    NumericVector out(nvox);
    out.attr("dim") = dims;
    double* op = REAL(out);

    // Uninitialised: each pass writes every element of its output before
    // anything reads it, so zero-filling these would be a wasted full-volume
    // sweep apiece. unique_ptr keeps the allocation exception-safe.
    std::unique_ptr<double[]> bufA_(new double[nvox]), bufB_(new double[nvox]);
    double* bufA = bufA_.get();
    double* bufB = bufB_.get();

    if (!normalize) {
        // Legacy path: a plain zero-padded convolution of the input, written
        // only at mask positions. Out-of-mask values are read, as before.
        const double* res = separable3(ap, bufA, bufB,
                                       d0, d1, d2, kx, ky, kz, window);
        if (full_mask) {
            std::memcpy(op, res, (size_t) nvox * sizeof(double));
            return out;
        }
        for (R_xlen_t i = 0; i < nm; i++) {
            const R_xlen_t v = (R_xlen_t) mask_idx[i] - 1;
            if (v < 0 || v >= nvox) stop("gaussian_blur_sep_cpp: 'mask_idx' out of bounds");
            op[v] = res[v];
        }
        return out;
    }

    // Mask-insulated path: out = blur(y) / blur(m).
    //
    // `full_mask` is the caller asserting that every voxel is in the mask, which
    // is what gaussian_blur(vol) with no mask means. Taking its word for it
    // skips materialising an nvox-long index vector in R and the indicator
    // scan here; when it is not set the indicator is built and counted as
    // before, and a full mask discovered that way takes the same shortcut.
    std::vector<char> in_mask;
    R_xlen_t n_distinct = nvox;
    if (!full_mask) {
        in_mask.assign((size_t) nvox, 0);
        n_distinct = 0;
        for (R_xlen_t i = 0; i < nm; i++) {
            const R_xlen_t v = (R_xlen_t) mask_idx[i] - 1;
            if (v < 0 || v >= nvox) stop("gaussian_blur_sep_cpp: 'mask_idx' out of bounds");
            if (!in_mask[v]) { in_mask[v] = 1; n_distinct++; }
        }
    }

    // When every voxel is in the mask -- which is what gaussian_blur(vol) with
    // no mask hands us, and the single most common call -- blur(m) is the
    // zero-padded convolution of the constant 1. That is separable in closed
    // form: the correction is a product of three per-axis factors, each just
    // the kernel weight that stays in bounds at that coordinate. So three
    // full-volume passes and an nvox buffer collapse to d0 + d1 + d2 additions,
    // and the numerator no longer needs the per-tap mask test either.
    if (n_distinct == nvox) {
        const std::vector<double> cx = edge_factor(d0, kx, window);
        const std::vector<double> cy = edge_factor(d1, ky, window);
        const std::vector<double> cz = edge_factor(d2, kz, window);

        const double* num = separable3(ap, bufA, bufB, d0, d1, d2, kx, ky, kz, window);
        const R_xlen_t slice = (R_xlen_t) d0 * d1;
        for (int z = 0; z < d2; z++) {
            for (int y = 0; y < d1; y++) {
                const double cyz = cy[y] * cz[z];
                const R_xlen_t base = (R_xlen_t) y * d0 + (R_xlen_t) z * slice;
                for (int x = 0; x < d0; x++) {
                    const double den = cx[x] * cyz;
                    op[base + x] = (den > 0.0) ? (num[base + x] / den) : 0.0;
                }
            }
        }
        return out;
    }

    // blur(m), landing in `denom`. The x pass reads the indicator directly, so
    // m is never materialised as doubles; ordering the ping-pong to finish in
    // `denom` means the result needs no copy out before bufA/bufB are reused.
    std::unique_ptr<double[]> denom_(new double[nvox]);
    double* denom = denom_.get();
    pass_x_masked(nullptr, in_mask.data(), bufB, d0, d1, d2, kx, window);
    pass_y(bufB, bufA, d0, d1, d2, ky, window);
    pass_z(bufA, denom, d0, d1, d2, kz, window);

    // blur(y), landing in bufA. Again the masking happens inside the x pass.
    pass_x_masked(ap, in_mask.data(), bufA, d0, d1, d2, kx, window);
    pass_y(bufA, bufB, d0, d1, d2, ky, window);
    pass_z(bufB, bufA, d0, d1, d2, kz, window);
    const double* num = bufA;

    for (R_xlen_t i = 0; i < nm; i++) {
        const R_xlen_t v = (R_xlen_t) mask_idx[i] - 1;
        // Unreachable in practice -- the centre voxel is in-mask by
        // construction, so its own kernel weight makes the denominator
        // positive -- but the dense path carries the same guard, and 0/0 here
        // would be a NaN rather than the dense path's 0.
        op[v] = (denom[v] > 0.0) ? (num[v] / denom[v]) : 0.0;
    }
    return out;
}

// ---------------------------------------------------------------------------
// 4-D driver
//
// Smoothing a run is the same kernel applied to each volume, and the R loop
// that used to do it paid for that structure three times over: it copied a
// volume out of the run and the result back in (two full passes over the data
// per volume, in R), it re-allocated the two nvox scratch buffers on every
// iteration, and -- with a partial mask -- it recomputed blur(m), which does
// not depend on the volume at all, once per volume.
//
// Everything invariant across volumes is hoisted here: the three axis kernels,
// the mask indicator, and either the closed-form edge correction (full mask) or
// the single blur(m) denominator (partial mask). What is left per volume is the
// three convolution passes and the divide, which is exactly the work, and it
// parallelises over volumes with no shared mutable state -- each worker owns
// its scratch buffers and writes only into its own volumes of the output.
//
// The per-volume arithmetic is unchanged and happens in the same order, so the
// result is bit-identical to smoothing the volumes one at a time. That is
// pinned in tests/testthat/test-gaussian-blur-separable.R.
// ---------------------------------------------------------------------------

namespace {

// Everything about the kernel that does not depend on which volume is being
// blurred. Held by reference in the worker; must outlive it.
struct BlurPlan4D {
    int d0, d1, d2, window;
    R_xlen_t nvox;
    std::vector<double> kx, ky, kz;
    bool normalize;
    bool whole_volume;              // every voxel in mask -> write everywhere
    std::vector<char> in_mask;      // empty when whole_volume
    std::vector<R_xlen_t> mpos;     // 0-based mask positions; empty when whole_volume
    std::vector<double> cx, cy, cz; // edge correction, whole_volume + normalize
    std::vector<double> denom;      // blur(m), !whole_volume + normalize
};

struct Blur4DWorker : public RcppParallel::Worker {
    const double* src;
    double* dst;
    const BlurPlan4D& P;

    Blur4DWorker(const double* src_, double* dst_, const BlurPlan4D& P_)
        : src(src_), dst(dst_), P(P_) {}

    void operator()(std::size_t begin, std::size_t end) {
        // Allocated once per chunk, not once per volume. Uninitialised for the
        // same reason as in the 3-D driver: every pass writes its whole output
        // before anything reads it.
        std::unique_ptr<double[]> a_(new double[P.nvox]), b_(new double[P.nvox]);
        double* a = a_.get();
        double* b = b_.get();
        const R_xlen_t slice = (R_xlen_t) P.d0 * P.d1;

        for (std::size_t t = begin; t < end; t++) {
            const double* ip = src + (R_xlen_t) t * P.nvox;
            double* op = dst + (R_xlen_t) t * P.nvox;

            if (!P.normalize) {
                const double* res = separable3(ip, a, b, P.d0, P.d1, P.d2,
                                               P.kx, P.ky, P.kz, P.window);
                if (P.whole_volume) {
                    std::memcpy(op, res, (size_t) P.nvox * sizeof(double));
                } else {
                    for (size_t i = 0; i < P.mpos.size(); i++) op[P.mpos[i]] = res[P.mpos[i]];
                }
                continue;
            }

            if (P.whole_volume) {
                const double* num = separable3(ip, a, b, P.d0, P.d1, P.d2,
                                               P.kx, P.ky, P.kz, P.window);
                for (int z = 0; z < P.d2; z++) {
                    for (int y = 0; y < P.d1; y++) {
                        const double cyz = P.cy[y] * P.cz[z];
                        const R_xlen_t base = (R_xlen_t) y * P.d0 + (R_xlen_t) z * slice;
                        for (int x = 0; x < P.d0; x++) {
                            const double den = P.cx[x] * cyz;
                            op[base + x] = (den > 0.0) ? (num[base + x] / den) : 0.0;
                        }
                    }
                }
                continue;
            }

            // Partial mask: blur(y) here, over the denominator computed once.
            pass_x_masked(ip, P.in_mask.data(), a, P.d0, P.d1, P.d2, P.kx, P.window);
            pass_y(a, b, P.d0, P.d1, P.d2, P.ky, P.window);
            pass_z(b, a, P.d0, P.d1, P.d2, P.kz, P.window);
            for (size_t i = 0; i < P.mpos.size(); i++) {
                const R_xlen_t v = P.mpos[i];
                op[v] = (P.denom[v] > 0.0) ? (a[v] / P.denom[v]) : 0.0;
            }
        }
    }
};

} // anonymous namespace

// Separable Gaussian blur of a 4-D image, volume by volume. Same arguments as
// gaussian_blur_sep_cpp() except that `arr` carries a 4-D dim attribute; the
// return carries the same one.
// [[Rcpp::export]]
NumericVector gaussian_blur_sep_4d_cpp(NumericVector arr, IntegerVector mask_idx,
                                       int window, double sigma,
                                       NumericVector spacing, bool normalize = true,
                                       bool full_mask = false) {
    SEXP dim_attr = arr.attr("dim");
    if (Rf_isNull(dim_attr) || Rf_xlength(dim_attr) != 4)
        stop("gaussian_blur_sep_4d_cpp: 'arr' must have a 4-D dim attribute");
    IntegerVector dims = dim_attr;
    const int d0 = dims[0], d1 = dims[1], d2 = dims[2], d3 = dims[3];
    if (d0 <= 0 || d1 <= 0 || d2 <= 0 || d3 <= 0)
        stop("gaussian_blur_sep_4d_cpp: 'dim' must be positive");
    if (window < 0) stop("gaussian_blur_sep_4d_cpp: 'window' must be non-negative");
    if (!(sigma > 0) || !R_finite(sigma)) stop("gaussian_blur_sep_4d_cpp: 'sigma' must be positive");
    if (spacing.size() < 3) stop("gaussian_blur_sep_4d_cpp: 'spacing' must have 3 elements");

    BlurPlan4D P;
    P.d0 = d0; P.d1 = d1; P.d2 = d2; P.window = window;
    P.nvox = (R_xlen_t) d0 * d1 * d2;
    P.normalize = normalize;
    if (P.nvox * (R_xlen_t) d3 > arr.size())
        stop("gaussian_blur_sep_4d_cpp: 'arr' is shorter than prod(dim)");

    P.kx = axis_kernel(window, sigma, spacing[0]);
    P.ky = axis_kernel(window, sigma, spacing[1]);
    P.kz = axis_kernel(window, sigma, spacing[2]);

    // Mask bookkeeping, and every bounds check, happens here rather than in the
    // worker: stop() unwinds through the R API and must not be called from a
    // parallel region.
    P.whole_volume = full_mask;
    if (!full_mask) {
        std::vector<char> m((size_t) P.nvox, 0);
        R_xlen_t n_distinct = 0;
        const R_xlen_t nm = mask_idx.size();
        for (R_xlen_t i = 0; i < nm; i++) {
            const R_xlen_t v = (R_xlen_t) mask_idx[i] - 1;
            if (v < 0 || v >= P.nvox) stop("gaussian_blur_sep_4d_cpp: 'mask_idx' out of bounds");
            if (!m[v]) { m[v] = 1; n_distinct++; }
        }
        if (n_distinct == P.nvox) {
            P.whole_volume = true;
        } else {
            P.in_mask.swap(m);
            // Distinct and ascending: the 3-D driver walks mask_idx as given and
            // may write a duplicate position twice, which is the same value.
            P.mpos.reserve((size_t) n_distinct);
            for (R_xlen_t v = 0; v < P.nvox; v++) if (P.in_mask[v]) P.mpos.push_back(v);
        }
    }

    if (normalize && P.whole_volume) {
        P.cx = edge_factor(d0, P.kx, window);
        P.cy = edge_factor(d1, P.ky, window);
        P.cz = edge_factor(d2, P.kz, window);
    } else if (normalize) {
        // blur(m) is the same for every volume, so it is computed once here
        // instead of d3 times inside the loop.
        P.denom.assign((size_t) P.nvox, 0.0);
        std::unique_ptr<double[]> a_(new double[P.nvox]), b_(new double[P.nvox]);
        pass_x_masked(nullptr, P.in_mask.data(), b_.get(), d0, d1, d2, P.kx, window);
        pass_y(b_.get(), a_.get(), d0, d1, d2, P.ky, window);
        pass_z(a_.get(), P.denom.data(), d0, d1, d2, P.kz, window);
    }

    NumericVector out((R_xlen_t) P.nvox * d3);
    out.attr("dim") = dims;

    Blur4DWorker worker(REAL(arr), REAL(out), P);
    RcppParallel::parallelFor(0, (std::size_t) d3, worker);
    return out;
}
