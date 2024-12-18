// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "indexFuns.h"
using namespace Rcpp;

// Calculate bilateral filter weights
// [[Rcpp::export]]
NumericMatrix bilateral_weights(int window, double spatial_sigma, double intensity_sigma, NumericVector spacing, double intensity_sd) {
  int sz = (2*window + 1);
  int total = sz * sz * sz;
  NumericMatrix out(total, 2);

  double spatial_denom = 2.0 * (spatial_sigma * spatial_sigma);
  double intensity_denom = 2.0 * (intensity_sigma * intensity_sigma * intensity_sd * intensity_sd);

  int ind = 0;
  for (int i = -window; i <= window; i++) {
    for (int j = -window; j <= window; j++) {
      for (int k = -window; k <= window; k++) {
        double dist2 = (i*i*spacing[0]*spacing[0]) + 
                       (j*j*spacing[1]*spacing[1]) + 
                       (k*k*spacing[2]*spacing[2]);
        double spatial_w = std::exp(-dist2 / spatial_denom);
        out(ind, 0) = spatial_w;
        out(ind, 1) = intensity_denom;
        ind++;
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector bilateral_filter_cpp(NumericVector arr, IntegerVector mask_idx, int window, double spatial_sigma, double intensity_sigma,
                                   NumericVector spacing) {
  IntegerVector dims = arr.attr("dim");
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];

  NumericVector out(arr.size());

  double intensity_sd = indexfuns::masked_sd(arr, mask_idx);
  NumericMatrix wts = bilateral_weights(window, spatial_sigma, intensity_sigma, spacing, intensity_sd);

  // Precompute intensity denominator (same for all)
  double intensity_denom = wts(0,1); 
  NumericVector spatial_wts = wts(_,0); 
  double *spatial_ptr = &spatial_wts[0];
  
  NumericMatrix cds = indexToGridCpp(mask_idx, dims);
  int slicedim = nx*ny;

  double *arr_ptr = &arr[0];
  double *out_ptr = &out[0];

  int N = spatial_wts.size();

  // Precompute neighborhood offsets once
  std::vector<int> offsets;
  offsets.reserve(N);

  {
    int ind = 0;
    for (int i = -window; i <= window; i++) {
      for (int j = -window; j <= window; j++) {
        for (int k = -window; k <= window; k++, ind++) {
          int offset = i + j*nx + k*slicedim;
          offsets.push_back(offset);
        }
      }
    }
  }

  for (int m = 0; m < mask_idx.size(); m++) {
    int x = cds(m,0)-1;
    int y = cds(m,1)-1;
    int z = cds(m,2)-1;
    int center_idx = x + y*nx + z*slicedim;
    double center_val = arr_ptr[center_idx];

    double val_sum = 0.0;
    double w_sum = 0.0;

    // Gather neighborhood
    for (int n = 0; n < N; n++) {
      int offset = offsets[n];
      // Compute neighbor coordinates to check bounds
      int k_vox = z + (offset/(nx*ny)); 
      int remain_k = offset % (nx*ny);
      int j_vox = y + (remain_k / nx);
      int i_vox = x + (remain_k % nx);

      double neigh_val = 0.0;
      if (i_vox >= 0 && i_vox < nx && 
          j_vox >= 0 && j_vox < ny &&
          k_vox >= 0 && k_vox < nz) {
        neigh_val = arr_ptr[center_idx + offset];
      }

      double diff = neigh_val - center_val;
      double diff2 = diff * diff;
      double w = spatial_ptr[n] * std::exp(-diff2 / intensity_denom);
      val_sum += w * neigh_val;
      w_sum += w;
    }

    out_ptr[mask_idx[m]-1] = val_sum / w_sum;
  }

  out.attr("dim") = dims;
  return out;
}




