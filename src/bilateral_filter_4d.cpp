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
#include <cmath>
#include "indexFuns.h"

using namespace Rcpp;
using namespace RcppParallel;

// Assume masked_sd is available:
// double masked_sd(NumericVector arr, IntegerVector mask_idx);

// [[Rcpp::export]]
NumericVector bilateral_filter_4d_cpp_par(
    NumericVector arr,
    IntegerVector mask_idx,
    int spatial_window,
    int temporal_window,
    double spatial_sigma,
    double intensity_sigma,
    double temporal_sigma,
    NumericVector spacing,
    double range_scale
) {
  // Extract dims from arr
  IntegerVector dims = arr.attr("dim");
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
  int nt = dims[3];

  int nvox = arr.size();
  NumericVector output(nvox);
  std::copy(arr.begin(), arr.end(), output.begin());

  int slicedim_xy = nx * ny;
  int slicedim_xyz = slicedim_xy * nz;
  double* arr_ptr = &arr[0];
  double* out_ptr = &output[0];

  // Compute intensity SD across all time points for masked voxels (robust to non-finite)
  long double sum = 0.0L;
  long double sumsq = 0.0L;
  size_t count = 0;
  for (int m = 0; m < mask_idx.size(); ++m) {
    int spatial_linear_idx = mask_idx[m] - 1; // zero-based
    for (int t = 0; t < nt; ++t) {
      double v = arr_ptr[spatial_linear_idx + t * slicedim_xyz];
      if (R_finite(v)) {
        sum += v;
        sumsq += static_cast<long double>(v) * static_cast<long double>(v);
        ++count;
      }
    }
  }
  double intensity_sd = range_scale;
  if (!std::isfinite(intensity_sd)) {
    intensity_sd = 0.0;
    if (count > 1) {
      long double mean = sum / static_cast<long double>(count);
      long double var = (sumsq / static_cast<long double>(count)) - mean * mean;
      intensity_sd = (var > 0.0L) ? std::sqrt(static_cast<double>(var)) : 0.0;
    }
  }

  // Precompute constants
  double spatial_var = 2.0 * spatial_sigma * spatial_sigma;
  double temporal_var = 2.0 * temporal_sigma * temporal_sigma;
  double intensity_var = 2.0 * intensity_sigma * intensity_sigma * intensity_sd * intensity_sd;
  const double min_intensity_var = 1e-12;
  if (!std::isfinite(intensity_var) || intensity_var < min_intensity_var) {
    intensity_var = min_intensity_var;
  }

  int sw = (spatial_window * 2) + 1;
  int tw = (temporal_window * 2) + 1;
  int total_elements = sw * sw * sw * tw;

  // Precompute spatial-temporal kernel
  std::vector<double> spatial_temporal_kernel(total_elements);
  {
    int idx = 0;
    for (int tt = -temporal_window; tt <= temporal_window; tt++) {
      double dt2 = (tt * spacing[3]) * (tt * spacing[3]);
      double w_temporal = std::exp(-dt2 / temporal_var);
      for (int ii = -spatial_window; ii <= spatial_window; ii++) {
        double dx2 = (ii * spacing[0]) * (ii * spacing[0]);
        for (int jj = -spatial_window; jj <= spatial_window; jj++) {
          double dy2 = (jj * spacing[1]) * (jj * spacing[1]);
          for (int kk = -spatial_window; kk <= spatial_window; kk++) {
            double dz2 = (kk * spacing[2]) * (kk * spacing[2]);
            double spatial_dist2 = dx2 + dy2 + dz2;
            double w_spatial = std::exp(-spatial_dist2 / spatial_var);
            spatial_temporal_kernel[idx] = w_spatial * w_temporal;
            idx++;
          }
        }
      }
    }
  }

  std::vector<unsigned char> in_mask;
  in_mask.assign(slicedim_xyz, 0);
  for (int m = 0; m < mask_idx.size(); ++m) {
    int idx = mask_idx[m] - 1;
    if (idx >= 0 && idx < slicedim_xyz) {
      in_mask[idx] = 1;
    }
  }

  // Precompute neighbor deltas
  std::vector<int> neighbor_dx(total_elements);
  std::vector<int> neighbor_dy(total_elements);
  std::vector<int> neighbor_dz(total_elements);
  std::vector<int> neighbor_dt(total_elements);
  {
    int idx = 0;
    for (int tt = -temporal_window; tt <= temporal_window; tt++) {
      for (int ii = -spatial_window; ii <= spatial_window; ii++) {
        for (int jj = -spatial_window; jj <= spatial_window; jj++) {
          for (int kk = -spatial_window; kk <= spatial_window; kk++) {
            neighbor_dx[idx] = ii;
            neighbor_dy[idx] = jj;
            neighbor_dz[idx] = kk;
            neighbor_dt[idx] = tt;
            idx++;
          }
        }
      }
    }
  }

  // A parallel worker
  struct BilateralFilter4DWorker : public Worker {
    const double* arr_ptr;
    double* out_ptr;
    const int nx, ny, nz, nt;
    const int slicedim_xy, slicedim_xyz;
    const std::vector<int>& neighbor_dx;
    const std::vector<int>& neighbor_dy;
    const std::vector<int>& neighbor_dz;
    const std::vector<int>& neighbor_dt;
    const std::vector<double>& spatial_temporal_kernel;
    const std::vector<unsigned char>& in_mask;
    const IntegerVector& mask_idx;
    double intensity_var;

    BilateralFilter4DWorker(
      const double* arr_ptr,
      double* out_ptr,
      int nx, int ny, int nz, int nt,
      int slicedim_xy, int slicedim_xyz,
      const std::vector<int>& neighbor_dx,
      const std::vector<int>& neighbor_dy,
      const std::vector<int>& neighbor_dz,
      const std::vector<int>& neighbor_dt,
      const std::vector<double>& spatial_temporal_kernel,
      const std::vector<unsigned char>& in_mask,
      const IntegerVector& mask_idx,
      double intensity_var
    ) : arr_ptr(arr_ptr), out_ptr(out_ptr), nx(nx), ny(ny), nz(nz), nt(nt),
        slicedim_xy(slicedim_xy), slicedim_xyz(slicedim_xyz),
        neighbor_dx(neighbor_dx), neighbor_dy(neighbor_dy), neighbor_dz(neighbor_dz),
        neighbor_dt(neighbor_dt), spatial_temporal_kernel(spatial_temporal_kernel),
        in_mask(in_mask), mask_idx(mask_idx), intensity_var(intensity_var) {}

    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t m = begin; m < end; m++) {
        int spatial_linear_idx = mask_idx[m] - 1; // zero-based
        int z = spatial_linear_idx / (nx*ny);
        int remain = spatial_linear_idx % (nx*ny);
        int y = remain / nx;
        int x = remain % nx;

        for (int t = 0; t < nt; t++) {
          int center_idx = x + y*nx + z*slicedim_xy + t*slicedim_xyz;
          double center_val = arr_ptr[center_idx];
          if (!R_finite(center_val)) {
            continue;
          }

          double val_sum = 0.0;
          double w_sum = 0.0;

          for (size_t n = 0; n < spatial_temporal_kernel.size(); n++) {
            int neighbor_t = t + neighbor_dt[n];
            int neighbor_z = z + neighbor_dz[n];
            int neighbor_y = y + neighbor_dy[n];
            int neighbor_x = x + neighbor_dx[n];

            // Check bounds
            if (neighbor_x < 0 || neighbor_x >= nx ||
                neighbor_y < 0 || neighbor_y >= ny ||
                neighbor_z < 0 || neighbor_z >= nz ||
                neighbor_t < 0 || neighbor_t >= nt) {
              continue;
            }

            int neighbor_spatial_idx = neighbor_x + neighbor_y*nx + neighbor_z*slicedim_xy;
            if (!in_mask[neighbor_spatial_idx]) {
              continue;
            }

            int neighbor_idx = neighbor_spatial_idx + neighbor_t*slicedim_xyz;
            double neighbor_val = arr_ptr[neighbor_idx];
            if (!R_finite(neighbor_val)) {
              continue;
            }
            double diff = (center_val - neighbor_val);
            double intensity_weight = std::exp(-(diff*diff) / intensity_var);

            double w = spatial_temporal_kernel[n] * intensity_weight;
            val_sum += w * neighbor_val;
            w_sum += w;
          }

          if (w_sum > 0.0) {
            out_ptr[center_idx] = val_sum / w_sum;
          } else {
            out_ptr[center_idx] = center_val;
          }
        }
      }
    }
  };

  BilateralFilter4DWorker worker(
    arr_ptr, out_ptr, nx, ny, nz, nt, slicedim_xy, slicedim_xyz,
    neighbor_dx, neighbor_dy, neighbor_dz, neighbor_dt,
    spatial_temporal_kernel, in_mask, mask_idx, intensity_var
  );

  // Parallel execution
  parallelFor(0, (size_t)mask_idx.size(), worker);

  output.attr("dim") = dims;
  return output;
}
