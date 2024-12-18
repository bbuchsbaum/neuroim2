#ifndef INDEXFUNS_H
#define INDEXFUNS_H

#include <Rcpp.h>

// Declare exported functions
Rcpp::NumericMatrix indexToGridCpp(Rcpp::IntegerVector idx, Rcpp::IntegerVector array_dim);
Rcpp::IntegerVector gridToIndex3DCpp(Rcpp::IntegerVector array_dim, Rcpp::NumericMatrix voxmat);
Rcpp::NumericVector gridToIndexCpp(Rcpp::IntegerVector array_dim, Rcpp::IntegerMatrix voxmat);
Rcpp::NumericVector exgridToIndex4DCpp(Rcpp::IntegerVector array_dim, 
                                      Rcpp::IntegerVector iind, 
                                      Rcpp::IntegerVector jind,
                                      Rcpp::IntegerVector kind, 
                                      Rcpp::IntegerVector mind);
Rcpp::NumericVector box_nbhd(Rcpp::NumericVector arr, 
                            Rcpp::IntegerVector dims, 
                            int x, int y, int z, 
                            int window, 
                            Rcpp::NumericVector out, 
                            int slicedim);
Rcpp::NumericVector gaussian_weights(int window, 
                                   double sigma, 
                                   Rcpp::NumericVector spacing);
Rcpp::NumericMatrix local_sphere(int vx, int vy, int vz, 
                                double radius, 
                                Rcpp::NumericVector spacing, 
                                Rcpp::IntegerVector dim);

namespace indexfuns {
    // Internal implementation functions
    Rcpp::NumericMatrix indexToGridCpp_impl(Rcpp::IntegerVector idx, Rcpp::IntegerVector array_dim);
    Rcpp::IntegerVector gridToIndex3DCpp_impl(Rcpp::IntegerVector array_dim, Rcpp::NumericMatrix voxmat);
    Rcpp::NumericVector gridToIndexCpp_impl(Rcpp::IntegerVector array_dim, Rcpp::IntegerMatrix voxmat);
    Rcpp::NumericVector exgridToIndex4DCpp_impl(Rcpp::IntegerVector array_dim, 
                                               Rcpp::IntegerVector iind, 
                                               Rcpp::IntegerVector jind,
                                               Rcpp::IntegerVector kind, 
                                               Rcpp::IntegerVector mind);
    
    Rcpp::NumericVector box_nbhd_impl(Rcpp::NumericVector arr, 
                                     Rcpp::IntegerVector dims, 
                                     int x, int y, int z, 
                                     int window, 
                                     Rcpp::NumericVector out, 
                                     int slicedim);
    
    Rcpp::NumericVector gaussian_weights_impl(int window, 
                                            double sigma, 
                                            Rcpp::NumericVector spacing);
    
    Rcpp::NumericVector gaussian_blur_cpp_impl(Rcpp::NumericVector arr, 
                                             Rcpp::IntegerVector mask_idx, 
                                             int window, 
                                             double sigma, 
                                             Rcpp::NumericVector spacing);
    
    Rcpp::NumericVector box_blur_impl(Rcpp::NumericVector arr, 
                                     Rcpp::IntegerVector mask_idx, 
                                     int window);
    
    Rcpp::NumericMatrix local_sphere_impl(int vx, int vy, int vz, 
                                         double radius, 
                                         Rcpp::NumericVector spacing, 
                                         Rcpp::IntegerVector dim);
    
    Rcpp::List local_spheres_impl(Rcpp::NumericMatrix centers, 
                                 double radius, 
                                 Rcpp::NumericVector spacing, 
                                 Rcpp::IntegerVector dim);

    // Helper functions
    inline double masked_sd(Rcpp::NumericVector arr, Rcpp::IntegerVector mask_idx) {
        int n = mask_idx.size();
        double sum = 0.0;
        double sum_squared = 0.0;

        for (int i = 0; i < n; i++) {
            double val = arr[mask_idx[i] - 1];
            sum += val;
            sum_squared += val * val;
        }

        double mean = sum / n;
        double variance = (sum_squared - (n * mean * mean)) / (n - 1);
        return sqrt(variance);
    }

    int coord3d_to_index(int x, int y, int z, int dx, int dy, int dz, int slicedim);
    int coord4d_to_index(int x, int y, int z, int t, int dim_x, int dim_y, int dim_z, int dim_t);
}

#endif