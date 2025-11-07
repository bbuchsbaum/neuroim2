// [[Rcpp::depends(Rcpp, RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

struct SpatialKernel {
  std::vector<int> offsets;
  std::vector<double> weights;
  int nx, ny, nz;
  int slice_xy;
  explicit SpatialKernel(int window,
                         double sigma_s,
                         const NumericVector& spacing,
                         int nx_,
                         int ny_,
                         int nz_)
      : nx(nx_), ny(ny_), nz(nz_), slice_xy(nx * ny) {
    int w = window;
    for (int dz = -w; dz <= w; ++dz) {
      for (int dy = -w; dy <= w; ++dy) {
        for (int dx = -w; dx <= w; ++dx) {
          int off = dx + dy * nx + dz * slice_xy;
          offsets.push_back(off);
          double dd = (dx * spacing[0]) * (dx * spacing[0]) +
                      (dy * spacing[1]) * (dy * spacing[1]) +
                      (dz * spacing[2]) * (dz * spacing[2]);
          double ws = std::exp(-(dd) / (2.0 * sigma_s * sigma_s));
          weights.push_back(ws);
        }
      }
    }
  }
};

inline double clip_r(double r) {
  if (!std::isfinite(r)) {
    return 0.0;
  }
  if (r > 0.999999) return 0.999999;
  if (r < -0.999999) return -0.999999;
  return r;
}

inline double corr_from_sums(long double sx,
                             long double sy,
                             long double sxx,
                             long double syy,
                             long double sxy,
                             long double n) {
  if (n < 3.0L) {
    return 0.0;
  }
  long double mx = sx / n;
  long double my = sy / n;
  long double vx = sxx / n - mx * mx;
  long double vy = syy / n - my * my;
  long double cov = sxy / n - mx * my;
  if (vx <= 0.0L || vy <= 0.0L) {
    return 0.0;
  }
  long double den = std::sqrt(vx * vy);
  if (den <= 0.0L) {
    return 0.0;
  }
  return clip_r(static_cast<double>(cov / den));
}

inline double corr_to_affinity(double r, int mode, double param) {
  if (mode == 0) {
    if (r <= 0.0) return 0.0;
    return std::pow(r, param);
  }
  if (mode == 1) {
    double u = 1.0 - r;
    return std::exp(-(u * u) / (2.0 * param * param));
  }
  double v = r - param;
  return (v > 0.0) ? v : 0.0;
}

inline double dot_self(const double* ptr, int len) {
  double acc = 0.0;
  for (int i = 0; i < len; ++i) {
    acc += ptr[i] * ptr[i];
  }
  return acc;
}

} // namespace

// ------------------------------------------------------------------
// Baseline CGB graph builder (no nuisance projection, unit weights)
// ------------------------------------------------------------------

// [[Rcpp::export]]
List build_cgb_graph_cpp(NumericVector arr,
                         IntegerVector mask_idx,
                         IntegerVector run_ends,
                         int window,
                         double spatial_sigma,
                         NumericVector spacing,
                         int corr_mode,
                         double corr_param,
                         int topk,
                         int leave_out_run,
                         NumericVector run_weights,
                         bool add_self) {
  IntegerVector dims = arr.attr("dim");
  if (dims.size() != 4) {
    stop("arr must be 4D [nx, ny, nz, nt]");
  }
  int nx = dims[0], ny = dims[1], nz = dims[2], nt = dims[3];
  int slice_xy = nx * ny;
  int slice_xyz = slice_xy * nz;
  int N = mask_idx.size();

  std::vector<uint8_t> in_mask(slice_xyz, 0);
  for (int i = 0; i < N; ++i) {
    int lin = mask_idx[i] - 1;
    if (lin >= 0 && lin < slice_xyz) {
      in_mask[lin] = 1;
    }
  }

  int K = run_ends.size();
  if (K < 1) {
    stop("run_ends must have length >= 1");
  }
  std::vector<int> r_start(K), r_end(K);
  int prev = 0;
  for (int k = 0; k < K; ++k) {
    r_start[k] = prev;
    r_end[k] = run_ends[k];
    prev = r_end[k];
  }
  if (prev != nt) {
    stop("run_ends do not sum to total timepoints");
  }

  std::vector<double> omega(K, 0.0);
  if (run_weights.size() == K) {
    for (int k = 0; k < K; ++k) {
      omega[k] = run_weights[k];
    }
  } else {
    for (int k = 0; k < K; ++k) {
      int nk = r_end[k] - r_start[k];
      omega[k] = std::max(1.0, static_cast<double>(nk) - 3.0);
    }
  }

  SpatialKernel kernel(window, spatial_sigma, spacing, nx, ny, nz);
  const double* data = arr.begin();

  auto count_row = [&](int idx) -> int {
    int lin = mask_idx[idx] - 1;
    int cz = lin / slice_xy;
    int rem = lin % slice_xy;
    int cy = rem / nx;
    int cx = rem % nx;

    std::vector<std::pair<double, int>> neigh;
    neigh.reserve(kernel.offsets.size());
    for (size_t m = 0; m < kernel.offsets.size(); ++m) {
      int off = kernel.offsets[m];
      int nz_ = cz + (off / slice_xy);
      int rem2 = off % slice_xy;
      int ny_ = cy + (rem2 / nx);
      int nx_ = cx + (rem2 % nx);
      if (nx_ < 0 || nx_ >= nx || ny_ < 0 || ny_ >= ny || nz_ < 0 || nz_ >= nz) continue;
      int n_lin = nx_ + ny_ * nx + nz_ * slice_xy;
      if (!in_mask[n_lin]) continue;
      if (add_self && n_lin == lin) continue;

      long double z_acc = 0.0L;
      long double w_acc = 0.0L;

      for (int k = 0; k < K; ++k) {
        if (leave_out_run >= 0 && k == leave_out_run) continue;
        int t0 = r_start[k];
        int t1 = r_end[k];
        long double sx = 0.0L, sy = 0.0L, sxx = 0.0L, syy = 0.0L, sxy = 0.0L, n = 0.0L;
        for (int t = t0; t < t1; ++t) {
          int ci = lin + t * slice_xyz;
          int ni = n_lin + t * slice_xyz;
          double xv = data[ci];
          double yv = data[ni];
          if (!R_finite(xv) || !R_finite(yv)) continue;
          sx += xv;
          sy += yv;
          sxx += static_cast<long double>(xv) * xv;
          syy += static_cast<long double>(yv) * yv;
          sxy += static_cast<long double>(xv) * yv;
          n += 1.0L;
        }
        double r = corr_from_sums(sx, sy, sxx, syy, sxy, n);
        double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
        z_acc += static_cast<long double>(omega[k] * z);
        w_acc += static_cast<long double>(omega[k]);
      }
      if (w_acc <= 0.0L) continue;
      double r_pool = std::tanh(static_cast<double>(z_acc / w_acc));
      double a_r = corr_to_affinity(r_pool, corr_mode, corr_param);
      if (a_r <= 0.0) continue;
      double w = a_r * kernel.weights[m];
      if (w > 0.0) {
        neigh.emplace_back(w, n_lin);
      }
    }

    if (neigh.empty() && !add_self) return 0;
    if (topk > 0 && static_cast<int>(neigh.size()) > topk) {
      std::nth_element(
        neigh.begin(),
        neigh.begin() + topk,
        neigh.end(),
        [](const std::pair<double,int>& lhs, const std::pair<double,int>& rhs) {
          return lhs.first > rhs.first;
        }
      );
      neigh.resize(topk);
    }
    double sumw = add_self ? 1e-12 : 0.0;
    for (const auto& p : neigh) sumw += p.first;
    if (sumw <= 0.0) return 0;
    return static_cast<int>(neigh.size()) + (add_self ? 1 : 0);
  };

  std::vector<int> nnz(N, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    nnz[i] = count_row(i);
  }

  std::vector<int> row_ptr(N + 1, 0);
  for (int i = 0; i < N; ++i) {
    row_ptr[i + 1] = row_ptr[i] + nnz[i];
  }
  int E = row_ptr[N];
  IntegerVector col_ind(E);
  NumericVector val(E);

  auto fill_row = [&](int idx, int write_pos) -> int {
    int lin = mask_idx[idx] - 1;
    int cz = lin / slice_xy;
    int rem = lin % slice_xy;
    int cy = rem / nx;
    int cx = rem % nx;

    std::vector<std::pair<double,int>> neigh;
    neigh.reserve(kernel.offsets.size());
    for (size_t m = 0; m < kernel.offsets.size(); ++m) {
      int off = kernel.offsets[m];
      int nz_ = cz + (off / slice_xy);
      int rem2 = off % slice_xy;
      int ny_ = cy + (rem2 / nx);
      int nx_ = cx + (rem2 % nx);
      if (nx_ < 0 || nx_ >= nx || ny_ < 0 || ny_ >= ny || nz_ < 0 || nz_ >= nz) continue;
      int n_lin = nx_ + ny_ * nx + nz_ * slice_xy;
      if (!in_mask[n_lin]) continue;
      if (add_self && n_lin == lin) continue;

      long double z_acc = 0.0L;
      long double w_acc = 0.0L;

      for (int k = 0; k < K; ++k) {
        if (leave_out_run >= 0 && k == leave_out_run) continue;
        int t0 = r_start[k];
        int t1 = r_end[k];
        long double sx = 0.0L, sy = 0.0L, sxx = 0.0L, syy = 0.0L, sxy = 0.0L, n = 0.0L;
        for (int t = t0; t < t1; ++t) {
          int ci = lin + t * slice_xyz;
          int ni = n_lin + t * slice_xyz;
          double xv = data[ci];
          double yv = data[ni];
          if (!R_finite(xv) || !R_finite(yv)) continue;
          sx += xv;
          sy += yv;
          sxx += static_cast<long double>(xv) * xv;
          syy += static_cast<long double>(yv) * yv;
          sxy += static_cast<long double>(xv) * yv;
          n += 1.0L;
        }
        double r = corr_from_sums(sx, sy, sxx, syy, sxy, n);
        double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
        z_acc += static_cast<long double>(omega[k] * z);
        w_acc += static_cast<long double>(omega[k]);
      }
      if (w_acc <= 0.0L) continue;
      double r_pool = std::tanh(static_cast<double>(z_acc / w_acc));
      double a_r = corr_to_affinity(r_pool, corr_mode, corr_param);
      if (a_r <= 0.0) continue;
      double w = a_r * kernel.weights[m];
      if (w > 0.0) {
        neigh.emplace_back(w, n_lin);
      }
    }

    if (neigh.empty() && !add_self) return 0;
    if (topk > 0 && static_cast<int>(neigh.size()) > topk) {
      std::nth_element(
        neigh.begin(),
        neigh.begin() + topk,
        neigh.end(),
        [](const std::pair<double,int>& lhs, const std::pair<double,int>& rhs) {
          return lhs.first > rhs.first;
        }
      );
      neigh.resize(topk);
    }
    // Merge self-edge (optional) with neighbors
    std::vector<std::pair<double,int>> cand;
    cand.reserve(neigh.size() + (add_self ? 1 : 0));
    if (add_self) cand.emplace_back(1e-12, lin);
    cand.insert(cand.end(), neigh.begin(), neigh.end());
    // Capacity promised in counting pass
    const int cap = row_ptr[idx + 1] - row_ptr[idx];
    if ((int)cand.size() > cap) {
      std::nth_element(
        cand.begin(),
        cand.begin() + cap,
        cand.end(),
        [](const std::pair<double,int>& L, const std::pair<double,int>& R){ return L.first > R.first; }
      );
      cand.resize(cap);
    }
    while ((int)cand.size() < cap) {
      cand.emplace_back(1e-12, lin);
    }
    double sumw = 0.0;
    for (const auto& p : cand) sumw += p.first;
    if (sumw <= 0.0) sumw = 1.0;
    const double inv = 1.0 / sumw;
    for (int k = 0; k < cap; ++k) {
      col_ind[write_pos + k] = cand[k].second;
      val[write_pos + k]     = cand[k].first * inv;
    }
    return cap;
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    fill_row(i, row_ptr[i]);
  }

  return List::create(
    _["row_ptr"] = wrap(row_ptr),
    _["col_ind"] = col_ind,
    _["val"] = val,
    _["dims3d"] = IntegerVector::create(nx, ny, nz),
    _["mask_idx"] = mask_idx
  );
}

// ------------------------------------------------------------------
// Builder with nuisance projection + time weights
// ------------------------------------------------------------------

// [[Rcpp::export]]
List build_cgb_graph_nuis_cpp(NumericVector arr,
                              IntegerVector mask_idx,
                              IntegerVector run_ends,
                              int window,
                              double spatial_sigma,
                              NumericVector spacing,
                              int corr_mode,
                              double corr_param,
                              int topk,
                              int leave_out_run,
                              NumericVector run_weights,
                              List Q_list,
                              List sqrtw_list,
                              bool add_self) {
  IntegerVector dims = arr.attr("dim");
  if (dims.size() != 4) {
    stop("arr must be 4D [nx, ny, nz, nt]");
  }
  int nx = dims[0], ny = dims[1], nz = dims[2], nt = dims[3];
  int slice_xy = nx * ny;
  int slice_xyz = slice_xy * nz;
  int N = mask_idx.size();

  std::vector<int> row_of(slice_xyz, -1);
  for (int i = 0; i < N; ++i) {
    int lin = mask_idx[i] - 1;
    if (lin >= 0 && lin < slice_xyz) {
      row_of[lin] = i;
    }
  }

  int K = run_ends.size();
  if (K < 1) stop("run_ends must have length >= 1");
  std::vector<int> r_start(K), r_end(K);
  int prev = 0;
  for (int k = 0; k < K; ++k) {
    r_start[k] = prev;
    r_end[k] = run_ends[k];
    prev = r_end[k];
  }
  if (prev != nt) {
    stop("run_ends do not sum to total timepoints");
  }

  std::vector<double> omega(K, 0.0);
  if (run_weights.size() == K) {
    for (int k = 0; k < K; ++k) omega[k] = run_weights[k];
  } else {
    for (int k = 0; k < K; ++k) {
      int nk = r_end[k] - r_start[k];
      omega[k] = std::max(1.0, static_cast<double>(nk) - 3.0);
    }
  }

  std::vector<bool> has_proj(K, false);
  std::vector<int> pcols(K, 0);
  for (int k = 0; k < K; ++k) {
    if (Q_list.size() == K && sqrtw_list.size() == K &&
        Q_list[k] != R_NilValue && sqrtw_list[k] != R_NilValue) {
      NumericMatrix Qk(Q_list[k]);
      pcols[k] = Qk.ncol();
      has_proj[k] = true;
    }
  }

  SpatialKernel kernel(window, spatial_sigma, spacing, nx, ny, nz);
  const double* data = arr.begin();

  size_t total_p = 0;
  for (int k = 0; k < K; ++k) {
    total_p += static_cast<size_t>(pcols[k]);
  }
  std::vector<size_t> p_offset(K + 1, 0);
  for (int k = 0; k < K; ++k) {
    p_offset[k + 1] = p_offset[k] + static_cast<size_t>(pcols[k]);
  }
  std::vector<double> sxx_w(static_cast<size_t>(K) * static_cast<size_t>(N), 0.0);
  std::vector<double> gx(static_cast<size_t>(N) * total_p, 0.0);

  std::vector<NumericMatrix> Qm(K);
  std::vector<NumericVector> Sw(K);
  std::vector<const double*> Q_ptr(K, nullptr);
  std::vector<const double*> Sw_ptr(K, nullptr);
  std::vector<int> Q_nrow(K, 0), Q_ncol(K, 0);
  for (int k = 0; k < K; ++k) {
    if (has_proj[k]) {
      Qm[k] = as<NumericMatrix>(Q_list[k]);
      Sw[k] = as<NumericVector>(sqrtw_list[k]);
      Q_ptr[k] = Qm[k].begin();
      Sw_ptr[k] = Sw[k].begin();
      Q_nrow[k] = Qm[k].nrow();
      Q_ncol[k] = Qm[k].ncol();
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int mi = 0; mi < N; ++mi) {
    int lin = mask_idx[mi] - 1;
    for (int k = 0; k < K; ++k) {
      int t0 = r_start[k];
      int t1 = r_end[k];
      long double sxx = 0.0L;
      if (has_proj[k]) {
        int pk = pcols[k];
        size_t off = p_offset[k];
        for (int p = 0; p < pk; ++p) {
          gx[static_cast<size_t>(mi) * total_p + off + p] = 0.0;
        }
        for (int t = 0; t < (t1 - t0); ++t) {
          int idx = lin + (t0 + t) * slice_xyz;
          double xv = data[idx];
          double sw = Sw_ptr[k][t];
          double xw = sw * xv;
          sxx += static_cast<long double>(xw) * xw;
          for (int p = 0; p < pk; ++p) {
            gx[static_cast<size_t>(mi) * total_p + off + p] +=
              Q_ptr[k][t + p * Q_nrow[k]] * xw;
          }
        }
      } else {
        for (int t = t0; t < t1; ++t) {
          int idx = lin + t * slice_xyz;
          double xv = data[idx];
          sxx += static_cast<long double>(xv) * xv;
        }
      }
      sxx_w[static_cast<size_t>(k) * static_cast<size_t>(N) + mi] = static_cast<double>(sxx);
    }
  }

  auto count_row = [&](int mi) -> int {
    int lin = mask_idx[mi] - 1;
    int cz = lin / slice_xy;
    int rem = lin % slice_xy;
    int cy = rem / nx;
    int cx = rem % nx;

    std::vector<std::pair<double,int>> neigh;
    neigh.reserve(kernel.offsets.size());

    for (size_t m = 0; m < kernel.offsets.size(); ++m) {
      int off = kernel.offsets[m];
      int nz_ = cz + (off / slice_xy);
      int rem2 = off % slice_xy;
      int ny_ = cy + (rem2 / nx);
      int nx_ = cx + (rem2 % nx);
      if (nx_ < 0 || nx_ >= nx || ny_ < 0 || ny_ >= ny || nz_ < 0 || nz_ >= nz) continue;
      int n_lin = nx_ + ny_ * nx + nz_ * slice_xy;
      int mj = (n_lin >= 0 && n_lin < slice_xyz) ? row_of[n_lin] : -1;
      if (mj < 0) continue;
      if (add_self && n_lin == lin) continue;

      long double z_acc = 0.0L;
      long double w_acc = 0.0L;

      for (int k = 0; k < K; ++k) {
        if (leave_out_run >= 0 && k == leave_out_run) continue;
        int t0 = r_start[k];
        int t1 = r_end[k];
        long double sxy = 0.0L;
        if (has_proj[k]) {
          for (int t = 0; t < (t1 - t0); ++t) {
            int ci = lin + (t0 + t) * slice_xyz;
            int ni = n_lin + (t0 + t) * slice_xyz;
            double w = Sw_ptr[k][t];
            sxy += static_cast<long double>(w * w) * data[ci] * data[ni];
          }
          size_t offp = p_offset[k];
          int pk = pcols[k];
          const double* gi = &gx[static_cast<size_t>(mi) * total_p + offp];
          const double* gj = &gx[static_cast<size_t>(mj) * total_p + offp];
          double dotg = 0.0;
          for (int p = 0; p < pk; ++p) {
            dotg += gi[p] * gj[p];
          }
          double sxx = sxx_w[static_cast<size_t>(k) * static_cast<size_t>(N) + mi];
          double syy = sxx_w[static_cast<size_t>(k) * static_cast<size_t>(N) + mj];
          double num = static_cast<double>(sxy) - dotg;
          double den_i = std::max(1e-12, sxx - dot_self(gi, pk));
          double den_j = std::max(1e-12, syy - dot_self(gj, pk));
          double den = std::sqrt(den_i * den_j);
          if (den <= 0.0) continue;
          double r = clip_r(num / den);
          double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
          z_acc += static_cast<long double>(omega[k] * z);
          w_acc += static_cast<long double>(omega[k]);
        } else {
          long double sxx = 0.0L, syy = 0.0L;
          for (int t = t0; t < t1; ++t) {
            int ci = lin + t * slice_xyz;
            int ni = n_lin + t * slice_xyz;
            double xi = data[ci];
            double xj = data[ni];
            sxy += static_cast<long double>(xi) * xj;
            sxx += static_cast<long double>(xi) * xi;
            syy += static_cast<long double>(xj) * xj;
          }
          double den = std::sqrt(static_cast<double>(sxx * syy));
          if (den <= 0.0) continue;
          double r = clip_r(static_cast<double>(sxy) / den);
          double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
          z_acc += static_cast<long double>(omega[k] * z);
          w_acc += static_cast<long double>(omega[k]);
        }
      }

      if (w_acc <= 0.0L) continue;
      double r_pool = std::tanh(static_cast<double>(z_acc / w_acc));
      double a_r = corr_to_affinity(r_pool, corr_mode, corr_param);
      if (a_r <= 0.0) continue;
      double w = a_r * kernel.weights[m];
      if (w > 0.0) neigh.emplace_back(w, n_lin);
    }

    if (neigh.empty() && !add_self) return 0;
    if (topk > 0 && static_cast<int>(neigh.size()) > topk) {
      std::nth_element(
        neigh.begin(),
        neigh.begin() + topk,
        neigh.end(),
        [](const std::pair<double,int>& lhs, const std::pair<double,int>& rhs) {
          return lhs.first > rhs.first;
        }
      );
      neigh.resize(topk);
    }
    double sumw = add_self ? 1e-12 : 0.0;
    for (const auto& p : neigh) sumw += p.first;
    if (sumw <= 0.0) return 0;
    return static_cast<int>(neigh.size()) + (add_self ? 1 : 0);
  };

  std::vector<int> nnz(N, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    nnz[i] = count_row(i);
  }

  std::vector<int> row_ptr(N + 1, 0);
  for (int i = 0; i < N; ++i) {
    row_ptr[i + 1] = row_ptr[i] + nnz[i];
  }
  int E = row_ptr[N];
  IntegerVector col_ind(E);
  NumericVector val(E);

  auto fill_row = [&](int mi, int write_pos) -> int {
    int lin = mask_idx[mi] - 1;
    int cz = lin / slice_xy;
    int rem = lin % slice_xy;
    int cy = rem / nx;
    int cx = rem % nx;

    std::vector<std::pair<double,int>> neigh;
    neigh.reserve(kernel.offsets.size());

    for (size_t m = 0; m < kernel.offsets.size(); ++m) {
      int off = kernel.offsets[m];
      int nz_ = cz + (off / slice_xy);
      int rem2 = off % slice_xy;
      int ny_ = cy + (rem2 / nx);
      int nx_ = cx + (rem2 % nx);
      if (nx_ < 0 || nx_ >= nx || ny_ < 0 || ny_ >= ny || nz_ < 0 || nz_ >= nz) continue;
      int n_lin = nx_ + ny_ * nx + nz_ * slice_xy;
      int mj = (n_lin >= 0 && n_lin < slice_xyz) ? row_of[n_lin] : -1;
      if (mj < 0) continue;

      long double z_acc = 0.0L;
      long double w_acc = 0.0L;

      for (int k = 0; k < K; ++k) {
        if (leave_out_run >= 0 && k == leave_out_run) continue;
        int t0 = r_start[k];
        int t1 = r_end[k];
        long double sxy = 0.0L;
        if (has_proj[k]) {
          for (int t = 0; t < (t1 - t0); ++t) {
            int ci = lin + (t0 + t) * slice_xyz;
            int ni = n_lin + (t0 + t) * slice_xyz;
            double w = Sw_ptr[k][t];
            sxy += static_cast<long double>(w * w) * data[ci] * data[ni];
          }
          size_t offp = p_offset[k];
          int pk = pcols[k];
          const double* gi = &gx[static_cast<size_t>(mi) * total_p + offp];
          const double* gj = &gx[static_cast<size_t>(mj) * total_p + offp];
          double dotg = 0.0;
          for (int p = 0; p < pk; ++p) {
            dotg += gi[p] * gj[p];
          }
          double sxx = sxx_w[static_cast<size_t>(k) * static_cast<size_t>(N) + mi];
          double syy = sxx_w[static_cast<size_t>(k) * static_cast<size_t>(N) + mj];
          double num = static_cast<double>(sxy) - dotg;
          double den_i = std::max(1e-12, sxx - dot_self(gi, pk));
          double den_j = std::max(1e-12, syy - dot_self(gj, pk));
          double den = std::sqrt(den_i * den_j);
          if (den <= 0.0) continue;
          double r = clip_r(num / den);
          double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
          z_acc += static_cast<long double>(omega[k] * z);
          w_acc += static_cast<long double>(omega[k]);
        } else {
          long double sxx = 0.0L, syy = 0.0L;
          for (int t = t0; t < t1; ++t) {
            int ci = lin + t * slice_xyz;
            int ni = n_lin + t * slice_xyz;
            double xi = data[ci];
            double xj = data[ni];
            sxy += static_cast<long double>(xi) * xj;
            sxx += static_cast<long double>(xi) * xi;
            syy += static_cast<long double>(xj) * xj;
          }
          double den = std::sqrt(static_cast<double>(sxx * syy));
          if (den <= 0.0) continue;
          double r = clip_r(static_cast<double>(sxy) / den);
          double z = 0.5 * std::log((1.0 + r) / (1.0 - r));
          z_acc += static_cast<long double>(omega[k] * z);
          w_acc += static_cast<long double>(omega[k]);
        }
      }

      if (w_acc <= 0.0L) continue;
      double r_pool = std::tanh(static_cast<double>(z_acc / w_acc));
      double a_r = corr_to_affinity(r_pool, corr_mode, corr_param);
      if (a_r <= 0.0) continue;
      double w = a_r * kernel.weights[m];
      if (w > 0.0) neigh.emplace_back(w, n_lin);
    }

    if (neigh.empty() && !add_self) return 0;
    if (topk > 0 && static_cast<int>(neigh.size()) > topk) {
      std::nth_element(
        neigh.begin(),
        neigh.begin() + topk,
        neigh.end(),
        [](const std::pair<double,int>& lhs, const std::pair<double,int>& rhs) {
          return lhs.first > rhs.first;
        }
      );
      neigh.resize(topk);
    }
    // Merge self-edge (optional) with neighbors
    std::vector<std::pair<double,int>> cand;
    cand.reserve(neigh.size() + (add_self ? 1 : 0));
    if (add_self) cand.emplace_back(1e-12, lin);
    cand.insert(cand.end(), neigh.begin(), neigh.end());
    // Capacity promised in counting pass
    const int cap = row_ptr[mi + 1] - row_ptr[mi];
    if ((int)cand.size() > cap) {
      std::nth_element(
        cand.begin(),
        cand.begin() + cap,
        cand.end(),
        [](const std::pair<double,int>& L, const std::pair<double,int>& R){ return L.first > R.first; }
      );
      cand.resize(cap);
    }
    while ((int)cand.size() < cap) {
      cand.emplace_back(1e-12, lin);
    }
    double sumw = 0.0;
    for (const auto& p : cand) sumw += p.first;
    if (sumw <= 0.0) sumw = 1.0;
    const double inv = 1.0 / sumw;
    for (int k = 0; k < cap; ++k) {
      col_ind[write_pos + k] = cand[k].second;
      val[write_pos + k]     = cand[k].first * inv;
    }
    return cap;
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    fill_row(i, row_ptr[i]);
  }

  return List::create(
    _["row_ptr"] = wrap(row_ptr),
    _["col_ind"] = col_ind,
    _["val"] = val,
    _["dims3d"] = IntegerVector::create(nx, ny, nz),
    _["mask_idx"] = mask_idx
  );
}

// ------------------------------------------------------------------
// Apply CSR graph to NeuroVol/NeuroVec
// ------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector apply_cgb_graph_cpp(NumericVector arr,
                                  IntegerVector row_ptr,
                                  IntegerVector col_ind,
                                  NumericVector val,
                                  IntegerVector mask_idx,
                                  int passes = 1,
                                  double lambda = 1.0) {
  IntegerVector dims = arr.attr("dim");
  if (dims.size() < 3 || dims.size() > 4) {
    stop("arr must be 3D or 4D");
  }
  int nx = dims[0], ny = dims[1], nz = dims[2];
  bool has_t = (dims.size() == 4);
  int nt = has_t ? dims[3] : 1;
  int slice_xy = nx * ny;
  int slice_xyz = slice_xy * nz;
  int N = mask_idx.size();

  NumericVector cur = clone(arr);
  NumericVector tmp(arr.size());

  for (int pass = 0; pass < passes; ++pass) {
    std::copy(cur.begin(), cur.end(), tmp.begin());
    double* x_in = cur.begin();
    double* x_out = tmp.begin();

    for (int t = 0; t < nt; ++t) {
      int t_off = t * slice_xyz;
      for (int i = 0; i < N; ++i) {
        int center_lin = mask_idx[i] - 1;
        int start = row_ptr[i];
        int end = row_ptr[i + 1];
        double acc = 0.0;
        for (int p = start; p < end; ++p) {
          int j_lin = col_ind[p];
          double w = val[p];
          acc += w * x_in[j_lin + t_off];
        }
        if (lambda >= 1.0) {
          x_out[center_lin + t_off] = acc;
        } else {
          x_out[center_lin + t_off] =
            (1.0 - lambda) * x_in[center_lin + t_off] + lambda * acc;
        }
      }
    }
    std::copy(tmp.begin(), tmp.end(), cur.begin());
  }

  cur.attr("dim") = dims;
  return cur;
}
