// [[Rcpp::depends(Rcpp)]]
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector downsample_4d_cpp(NumericVector arr, 
                                IntegerVector new_dims,
                                IntegerVector old_dims) {
  
  // Validate dimension vectors
  if (old_dims.size() != 4 || new_dims.size() != 4) {
    stop("Both old_dims and new_dims must have exactly 4 elements");
  }
  
  // old_dims and new_dims are both 4D (x, y, z, t)
  int old_nx = old_dims[0];
  int old_ny = old_dims[1];
  int old_nz = old_dims[2];
  int nt = old_dims[3]; // time dimension stays the same
  
  int new_nx = new_dims[0];
  int new_ny = new_dims[1];
  int new_nz = new_dims[2];
  
  // Validate dimensions are positive
  if (old_nx <= 0 || old_ny <= 0 || old_nz <= 0 || nt <= 0) {
    stop("Old dimensions must all be positive");
  }
  if (new_nx <= 0 || new_ny <= 0 || new_nz <= 0) {
    stop("New dimensions must all be positive");
  }
  
  // Validate array size matches dimensions
  R_xlen_t expected_size = (R_xlen_t)old_nx * old_ny * old_nz * nt;
  if (arr.size() != expected_size) {
    stop("Input array size (" + std::to_string(arr.size()) + 
         ") does not match old_dims (" + std::to_string(expected_size) + ")");
  }
  
  // Calculate downsampling factors
  double factor_x = (double)old_nx / new_nx;
  double factor_y = (double)old_ny / new_ny;
  double factor_z = (double)old_nz / new_nz;
  
  // Create output array
  R_xlen_t out_size = (R_xlen_t)new_nx * new_ny * new_nz * nt;
  NumericVector output(out_size);
  
  // Precompute strides for old array (use R_xlen_t to prevent overflow)
  R_xlen_t old_stride_y = old_nx;
  R_xlen_t old_stride_z = (R_xlen_t)old_nx * old_ny;
  R_xlen_t old_stride_t = (R_xlen_t)old_nx * old_ny * old_nz;
  
  // Precompute strides for new array
  R_xlen_t new_stride_y = new_nx;
  R_xlen_t new_stride_z = (R_xlen_t)new_nx * new_ny;
  R_xlen_t new_stride_t = (R_xlen_t)new_nx * new_ny * new_nz;
  
  // For each time point
  for (int t = 0; t < nt; t++) {
    R_xlen_t old_t_offset = (R_xlen_t)t * old_stride_t;
    R_xlen_t new_t_offset = (R_xlen_t)t * new_stride_t;
    
    // For each output voxel
    for (int new_z = 0; new_z < new_nz; new_z++) {
      for (int new_y = 0; new_y < new_ny; new_y++) {
        for (int new_x = 0; new_x < new_nx; new_x++) {
          
          // Calculate the range of input voxels to average
          double old_x_start = new_x * factor_x;
          double old_x_end = (new_x + 1) * factor_x;
          double old_y_start = new_y * factor_y;
          double old_y_end = (new_y + 1) * factor_y;
          double old_z_start = new_z * factor_z;
          double old_z_end = (new_z + 1) * factor_z;
          
          // Convert to integer bounds
          int x_min = (int)std::floor(old_x_start);
          int x_max = std::min((int)std::ceil(old_x_end), old_nx);
          int y_min = (int)std::floor(old_y_start);
          int y_max = std::min((int)std::ceil(old_y_end), old_ny);
          int z_min = (int)std::floor(old_z_start);
          int z_max = std::min((int)std::ceil(old_z_end), old_nz);
          
          // Accumulate values in the box
          double sum = 0.0;
          double weight_sum = 0.0;
          
          for (int z = z_min; z < z_max; z++) {
            for (int y = y_min; y < y_max; y++) {
              for (int x = x_min; x < x_max; x++) {
                // Calculate partial volume weights for boundary voxels
                double weight_x = 1.0;
                double weight_y = 1.0;
                double weight_z = 1.0;
                
                // Adjust weights for partial voxels
                if (x == x_min && old_x_start > x_min) {
                  weight_x = 1.0 - (old_x_start - x_min);
                } else if (x == x_max - 1 && old_x_end < x_max) {
                  weight_x = old_x_end - x;
                }
                
                if (y == y_min && old_y_start > y_min) {
                  weight_y = 1.0 - (old_y_start - y_min);
                } else if (y == y_max - 1 && old_y_end < y_max) {
                  weight_y = old_y_end - y;
                }
                
                if (z == z_min && old_z_start > z_min) {
                  weight_z = 1.0 - (old_z_start - z_min);
                } else if (z == z_max - 1 && old_z_end < z_max) {
                  weight_z = old_z_end - z;
                }
                
                double weight = weight_x * weight_y * weight_z;
                
                R_xlen_t old_idx = old_t_offset + x + y * old_stride_y + z * old_stride_z;
                double val = arr[old_idx];
                
                // Skip NaN and Inf values
                if (R_finite(val)) {
                  sum += val * weight;
                  weight_sum += weight;
                }
              }
            }
          }
          
          // Store averaged value
          R_xlen_t new_idx = new_t_offset + new_x + new_y * new_stride_y + new_z * new_stride_z;
          output[new_idx] = (weight_sum > 0) ? (sum / weight_sum) : NA_REAL;
        }
      }
    }
  }
  
  // Set dimensions attribute
  output.attr("dim") = new_dims;
  
  return output;
}
