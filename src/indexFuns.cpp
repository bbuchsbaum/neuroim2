// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "indexFuns.h"
using namespace Rcpp;

// Exported functions (outside namespace)
// [[Rcpp::export]]
NumericMatrix indexToGridCpp(IntegerVector idx, IntegerVector array_dim) {
    return indexfuns::indexToGridCpp_impl(idx, array_dim);
}

// [[Rcpp::export]]
IntegerVector gridToIndex3DCpp(IntegerVector array_dim, NumericMatrix voxmat) {
    return indexfuns::gridToIndex3DCpp_impl(array_dim, voxmat);
}

// [[Rcpp::export]]
NumericVector gridToIndexCpp(IntegerVector array_dim, IntegerMatrix voxmat) {
    return indexfuns::gridToIndexCpp_impl(array_dim, voxmat);
}

// [[Rcpp::export]]
NumericVector exgridToIndex4DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind,
                                IntegerVector kind, IntegerVector mind) {
    return indexfuns::exgridToIndex4DCpp_impl(array_dim, iind, jind, kind, mind);
}

// [[Rcpp::export]]
NumericVector box_nbhd(NumericVector arr, IntegerVector dims, int x, int y, int z, 
                      int window, NumericVector out, int slicedim) {
    return indexfuns::box_nbhd_impl(arr, dims, x, y, z, window, out, slicedim);
}

// [[Rcpp::export]]
NumericVector gaussian_weights(int window, double sigma, NumericVector spacing) {
    return indexfuns::gaussian_weights_impl(window, sigma, spacing);
}

// [[Rcpp::export]]
NumericVector gaussian_blur_cpp(NumericVector arr, IntegerVector mask_idx, int window, 
                              double sigma, NumericVector spacing) {
    return indexfuns::gaussian_blur_cpp_impl(arr, mask_idx, window, sigma, spacing);
}

// [[Rcpp::export]]
NumericVector box_blur(NumericVector arr, IntegerVector mask_idx, int window) {
    return indexfuns::box_blur_impl(arr, mask_idx, window);
}

// [[Rcpp::export]]
NumericMatrix local_sphere(int vx, int vy, int vz, double radius, 
                         NumericVector spacing, IntegerVector dim) {
    return indexfuns::local_sphere_impl(vx, vy, vz, radius, spacing, dim);
}

// [[Rcpp::export]]
List local_spheres(NumericMatrix centers, double radius, 
                  NumericVector spacing, IntegerVector dim) {
    return indexfuns::local_spheres_impl(centers, radius, spacing, dim);
}

// Implementation namespace
namespace indexfuns {

    NumericMatrix indexToGridCpp_impl(IntegerVector idx, IntegerVector array_dim) {
        int rank = array_dim.size();
        int N = idx.size();
        NumericMatrix omat(idx.size(), array_dim.size());

        for(int i = 0; i < N; i++) {
            int wh1 = idx(i)-1;
            int tmp = 1 + wh1 % array_dim(0);
            IntegerVector wh = IntegerVector(rank, tmp);
            if (rank >= 2) {
                int denom = 1;
                for (int j = 1; j < rank; j++) {
                    denom = denom * array_dim(j-1);
                    int nextd1 = (int)wh1/denom;
                    wh(j) = 1 + nextd1 % array_dim(j);
                }
            }
            omat.row(i) = wh;
        }
        return omat;
    }

    IntegerVector gridToIndex3DCpp_impl(IntegerVector array_dim, NumericMatrix voxmat) {
        int slicedim = array_dim[0]*array_dim[1];
        IntegerVector out = IntegerVector(voxmat.nrow());

        for (int i=0; i < voxmat.nrow(); i++) {
            out[i] = (int)((slicedim * (voxmat(i,2) -1)) + ((voxmat(i,1)-1) * array_dim(0)) + voxmat(i,0));
        }
        return out;
    }

    NumericVector gridToIndexCpp_impl(IntegerVector array_dim, IntegerMatrix voxmat) {
        int ndim = array_dim.size();
        std::vector<int64_t> multipliers(ndim);
        multipliers[0] = 1;
        for (int i = 1; i < ndim; ++i) {
            multipliers[i] = multipliers[i - 1] * static_cast<int64_t>(array_dim[i - 1]);
        }

        int64_t nvox = voxmat.nrow();
        NumericVector out(nvox);

        for (int64_t i = 0; i < nvox; ++i) {
            int64_t ind = 0;
            for (int j = 0; j < ndim; ++j) {
                ind += (static_cast<int64_t>(voxmat(i, j)) - 1) * multipliers[j];
            }
            out[i] = static_cast<double>(ind + 1);
        }
        return out;
    }

    NumericVector exgridToIndex4DCpp_impl(IntegerVector array_dim, IntegerVector iind,
                                        IntegerVector jind, IntegerVector kind, IntegerVector mind) {
        int64_t X = array_dim[0];
        int64_t Y = array_dim[1];
        int64_t Z = array_dim[2];

        int64_t m1 = 1;
        int64_t m2 = X;
        int64_t m3 = X * Y;
        int64_t m4 = X * Y * Z;

        int64_t nels = static_cast<int64_t>(iind.size()) * jind.size() * kind.size() * mind.size();
        NumericVector out(nels);

        int64_t count = 0;
        for (int64_t m = 0; m < mind.size(); ++m) {
            int64_t m_offset = (mind[m] - 1) * m4;
            for (int64_t k = 0; k < kind.size(); ++k) {
                int64_t km_offset = m_offset + (kind[k] - 1) * m3;
                for (int64_t j = 0; j < jind.size(); ++j) {
                    int64_t jkm_offset = km_offset + (jind[j] - 1) * m2;
                    for (int64_t i = 0; i < iind.size(); ++i) {
                        int64_t ind = jkm_offset + (iind[i] - 1) * m1;
                        out[count++] = static_cast<double>(ind + 1);
                    }
                }
            }
        }
        return out;
    }

    NumericVector box_nbhd_impl(NumericVector arr, IntegerVector dims, int x, int y, int z, 
                               int window, NumericVector out, int slicedim) {
        int ind = 0;
        for (int k = z-window; k <= z+window; k++) {
            for (int j = y-window; j <= y+window; j++) {
                for (int i = x-window; i <= x+window; i++) {
                    if (i >= 0 && i < dims[0] && j >= 0 && j < dims[1] && k >= 0 && k < dims[2]) {
                        out[ind] = arr[i + j*dims[0] + k*slicedim];
                    } else {
                        out[ind] = 0;
                    }
                    ind++;
                }
            }
        }
        return out;
    }

    NumericVector gaussian_weights_impl(int window, double sigma, NumericVector spacing) {
        int sz = (2*window + 1);
        NumericVector out(sz*sz*sz);
        double denom = 2.0 * sigma * sigma;
        int ind = 0;

        for (int k = -window; k <= window; k++) {
            for (int j = -window; j <= window; j++) {
                for (int i = -window; i <= window; i++) {
                    out[ind] = std::exp(-std::pow(i * spacing[0],2)/denom)*std::sqrt(2 * M_PI * sigma) *
                              std::exp(-std::pow(j * spacing[1],2)/denom)*std::sqrt(2 * M_PI * sigma) *
                              std::exp(-std::pow(k * spacing[2],2)/denom)*std::sqrt(2 * M_PI * sigma);
                    ind++;
                }
            }
        }

        double tot = std::accumulate(out.begin(), out.end(), 0.0);
        return out/tot;
    }

    NumericVector gaussian_blur_cpp_impl(NumericVector arr, IntegerVector mask_idx, int window, 
                                       double sigma, NumericVector spacing) {
        IntegerVector dims = arr.attr("dim");
        NumericVector out = NumericVector(arr.length());
        NumericVector local = NumericVector(pow((window*2)+1, 3));

        NumericVector wts = gaussian_weights(window, sigma, spacing);
        NumericMatrix cds = indexToGridCpp(mask_idx, dims);
        int slicedim = dims[0]*dims[1];

        for (int i = 0; i < mask_idx.length(); i++) {
            NumericVector ret = box_nbhd(arr, dims, cds(i,0)-1, cds(i,1)-1, cds(i,2)-1, 
                                       window, local, slicedim);
            out[mask_idx[i]-1] = sum(ret*wts);
        }

        out.attr("dim") = dims;
        return out;
    }

    NumericVector box_blur_impl(NumericVector arr, IntegerVector mask_idx, int window) {
        IntegerVector dims = arr.attr("dim");
        NumericVector out = NumericVector(arr.length());
        NumericVector local = NumericVector(pow((window*2)+1, 3));

        NumericMatrix cds = indexToGridCpp(mask_idx, dims);
        int slicedim = dims[0]*dims[1];

        for (int i = 0; i < mask_idx.length(); i++) {
            NumericVector ret = box_nbhd(arr, dims, cds(i,0)-1, cds(i,1)-1, cds(i,2)-1, 
                                       window, local, slicedim);
            out[mask_idx[i]-1] = mean(ret);
        }

        out.attr("dim") = dims;
        return out;
    }

    NumericMatrix local_sphere_impl(int vx, int vy, int vz, double radius, 
                                  NumericVector spacing, IntegerVector dim) {
        std::vector<double> xs, ys, zs;
        int window = ceil(radius/min(spacing));
        
        for(int i = -window; i <= window; i++) {
            for(int j = -window; j <= window; j++) {
                for(int k = -window; k <= window; k++) {
                    double dist = sqrt(pow(i*spacing[0],2) + 
                                     pow(j*spacing[1],2) + 
                                     pow(k*spacing[2],2));
                    if(dist <= radius) {
                        int newx = vx + i;
                        int newy = vy + j;
                        int newz = vz + k;
                        if(newx >= 0 && newx < dim[0] &&
                           newy >= 0 && newy < dim[1] &&
                           newz >= 0 && newz < dim[2]) {
                            xs.push_back(newx);
                            ys.push_back(newy);
                            zs.push_back(newz);
                        }
                    }
                }
            }
        }
        
        NumericMatrix out(xs.size(), 3);
        for(size_t i = 0; i < xs.size(); i++) {
            out(i,0) = xs[i];
            out(i,1) = ys[i];
            out(i,2) = zs[i];
        }
        return out;
    }

    List local_spheres_impl(NumericMatrix centers, double radius, 
                           NumericVector spacing, IntegerVector dim) {
        List out(centers.nrow());
        for(int i = 0; i < centers.nrow(); i++) {
            out[i] = local_sphere(centers(i,0)-1, centers(i,1)-1, centers(i,2)-1,
                                radius, spacing, dim);
        }
        return out;
    }

    int coord3d_to_index(int x, int y, int z, int dx, int dy, int dz, int slicedim) {
        return (slicedim * z) + (y * dx) + x;
    }

    int coord4d_to_index(int x, int y, int z, int t, int dim_x, int dim_y, int dim_z, int dim_t) {
        if (x < 0 || x >= dim_x || y < 0 || y >= dim_y || z < 0 || z >= dim_z || t < 0 || t >= dim_t) {
            return -1;
        }
        return x + dim_x * (y + dim_y * (z + dim_z * t));
    }

} // end namespace indexfuns

