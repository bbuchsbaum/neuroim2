#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix indexToGridCpp(IntegerVector idx, IntegerVector array_dim) {
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



// NumericMatrix nn3d(IntegerVector pt, IntegerVector spacing, IntegerVector dim, double radius) {
//   int i_rad = floor(radius/spacing[0]);
//   int j_rad = floor(radius/spacing[1]);
//   int k_rad = floor(radius/spacing[2]);
//
//   NumericMatrix out((i_rad*2+1)*(j_rad*2+1)*(k_rad*2+1), 3);
//   NumericVector cd = NumericVector::create(pt[0]*spacing[0], pt[1]*spacing[1], pt[2]*spacing[2]);
//
//   for (int i = -i_rad; i <= i_rad; i++) {
//     for (int j = -j_rad; j <= j_rad; j++) {
//       for (int k = -k_rad; k <= k_rad; k++) {
//         double x = i*spacing[0];
//         double y = j*spacing[1];
//         double z = k*spacing[2];
//
//         double d = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
//
//       }
//     }
//   }
//
// }

// [[Rcpp::export]]
NumericMatrix local_sphere(int vx, int vy, int vz, double radius, NumericVector spacing, IntegerVector dim) {

  if (radius <= 0) {
    Rcpp::stop("radius must be greater than 0.");
  }

  if (spacing.size() != 3) {
    Rcpp::stop("spacing must be a vector of length 3.");
  }

  if (dim.size() != 3) {
    Rcpp::stop("dim must be an integer vector of length 3.");
  }

  int i_rad = round(radius/spacing[0]);
  int j_rad = round(radius/spacing[1]);
  int k_rad = round(radius/spacing[2]);

  //std::list<std::vector<double>> cds;
  //List cds((i_rad*2)*(j_rad*2)*(k_rad*2));

  NumericMatrix cds((i_rad*2+1)*(j_rad*2+1)*(k_rad*2+1), 3);

  int count = 0;
  for (int i = vx - i_rad; i <= (vx + i_rad); i++) {
      if (i < 1 || i > dim[0]) {
        continue;
      }
      //Rcout << "i: " << i << std::endl;
      for (int j = vy - j_rad; j <= (vy + j_rad); j++) {
        if (j < 1 || j > dim[1]) {
          continue;
        }
        //Rcout << "j: " << j << std::endl;
        for (int k = vz -k_rad; k <= (vz + k_rad); k++) {
          if (k < 1 || k > dim[2]) {
            continue;
          }
          //Rcout << "k: " << k << std::endl;

          double xd = (i-vx)*spacing[0];
          double yd = (j-vy)*spacing[1];
          double zd = (k-vz)*spacing[2];

          double d = sqrt(pow(xd,2) + pow(yd,2) + pow(zd,2));
          //Rcout << "d: " << d << std::endl;
          if (d < radius) {
            //cds[count] = IntegerVector::create(i,j,k);
            //count++;
            cds(count,0) = i;
            cds(count,1) = j;
            cds(count,2) = k;
            count++;
          }

        }
      }
  }
  return cds(Rcpp::Range(0, count-1), Rcpp::_);

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

// [[Rcpp::export]]
NumericVector box_nbhd(NumericVector arr, IntegerVector dims, int x, int y, int z, int window, NumericVector out, int slicedim) {
  int count = 0;
  //NumericVector out = NumericVector( pow((window*2)+1, 3));
  for (int i = x - window; i <= x + window; i++) {
    for (int j = y - window; j <= y + window; j++) {
      for (int k = z - window; k <= z + window; k++) {
        int ind = coord3d_to_index(i,j,k, dims[0], dims[1], dims[2], slicedim);
        if (ind < 0 || ind >= arr.length()) {
          out[count] = 0;
        } else {
          out[count] = arr[ind];
        }
        count = count +1;
      }
    }
  }

  return out;

}

NumericVector box_nbhd_4d(NumericVector arr, IntegerVector dims, int x, int y, int z, int t, int spatial_window, int temporal_window, NumericVector out) {
  int count = 0;

  for (int i = x - spatial_window; i <= x + spatial_window; i++) {
    for (int j = y - spatial_window; j <= y + spatial_window; j++) {
      for (int k = z - spatial_window; k <= z + spatial_window; k++) {
        for (int l = t - temporal_window; l <= t + temporal_window; l++) {
          int ind = coord4d_to_index(i, j, k, l, dims[0], dims[1], dims[2], dims[3]);
          if (ind < 0 || ind >= arr.length()) {
            out[count] = 0;
          } else {
            out[count] = arr[ind];
          }
          count = count + 1;
        }
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector box_blur(NumericVector arr, IntegerVector mask_idx, int window) {
  IntegerVector dims = arr.attr("dim");
  NumericVector out = NumericVector(arr.length());
  NumericVector local = NumericVector( pow((window*2)+1, 3));

  NumericMatrix cds = indexToGridCpp(mask_idx, dims);
  int slicedim = dims[0]*dims[1];

  for (int i = 0; i < mask_idx.length(); i++) {
    NumericVector ret = box_nbhd(arr, dims, cds(i,0)-1, cds(i,1)-1, cds(i,2)-1, window, local, slicedim);
    out[mask_idx[i]-1] = mean(ret);
  }

  out.attr("dim") = dims;
  return out;
}

// H(x) = exp(-x2/ (2s2)) / sqrt(2* pi*s2)
// I(y) = exp(-y2/ (2t2)) / sqrt(2* pi*t2)

// [[Rcpp::export]]
NumericVector gaussian_weights(int window, double sigma, NumericVector spacing) {
  int count = 0;
  NumericVector out = NumericVector( pow((window*2)+1, 3));
  double denom = 2 * pow(sigma,2);

  int ind = 0;
  for (int i = -window; i <= window; i++) {
    for (int j = -window; j <= window; j++) {
      for (int k = -window; k <= window; k++) {
          out[ind] = std::exp(-std::pow(i * spacing[0],2)/denom)*std::sqrt(2 * M_PI * sigma) *
            std::exp(-std::pow(j * spacing[1],2)/denom)*std::sqrt(2 * M_PI * sigma) *
            std::exp(-std::pow(k * spacing[2],2)/denom)*std::sqrt(2 * M_PI * sigma);
          ind = ind +1;
      }
    }
  }

  double tot = std::accumulate(out.begin(), out.end(), 0.0);
  return out/tot;


}


// [[Rcpp::export]]
NumericVector gaussian_blur_cpp(NumericVector arr, IntegerVector mask_idx, int window, double sigma, NumericVector spacing) {
  IntegerVector dims = arr.attr("dim");
  NumericVector out = NumericVector(arr.length());
  NumericVector local = NumericVector( pow((window*2)+1, 3));

  NumericVector wts = gaussian_weights(window, sigma, spacing);

  //Rcout << "weights are" << std::endl << wts << std::endl;

  NumericMatrix cds = indexToGridCpp(mask_idx, dims);
  int slicedim = dims[0]*dims[1];

  for (int i = 0; i < mask_idx.length(); i++) {
    NumericVector ret = box_nbhd(arr, dims, cds(i,0)-1, cds(i,1)-1, cds(i,2)-1, window, local, slicedim);
    out[mask_idx[i]-1] = sum(ret*wts);
  }

  out.attr("dim") = dims;
  return out;
}

// Compute standard deviation of intensity values within the mask
double masked_sd(NumericVector arr, IntegerVector mask_idx) {
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


// Calculate bilateral filter weights
NumericMatrix bilateral_weights(int window, double spatial_sigma, double intensity_sigma, NumericVector spacing, double intensity_sd) {
  NumericMatrix out = NumericMatrix(pow((window*2)+1, 3), 2);
  double spatial_denom = 2 * pow(spatial_sigma, 2);
  double intensity_denom = 2 * pow(intensity_sigma * intensity_sd, 2);

  int ind = 0;
  for (int i = -window; i <= window; i++) {
    for (int j = -window; j <= window; j++) {
      for (int k = -window; k <= window; k++) {
        out(ind, 0) = std::exp(-std::pow(i * spacing[0],2)/spatial_denom) * std::exp(-std::pow(j * spacing[1],2)/spatial_denom) * std::exp(-std::pow(k * spacing[2],2)/spatial_denom);
        out(ind, 1) = intensity_denom;
        ind = ind + 1;
      }
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericVector bilateral_filter_cpp(NumericVector arr, IntegerVector mask_idx, int window, double spatial_sigma, double intensity_sigma,
                                   NumericVector spacing) {
  IntegerVector dims = arr.attr("dim");
  NumericVector out = NumericVector(arr.length());
  NumericVector local = NumericVector( pow((window*2)+1, 3));

  // Find the standard deviation of the intensity values within the mask
  double intensity_sd = masked_sd(arr, mask_idx);

  // Calculate bilateral weights
  NumericMatrix wts = bilateral_weights(window, spatial_sigma, intensity_sigma, spacing, intensity_sd);

  NumericMatrix cds = indexToGridCpp(mask_idx, dims);
  int slicedim = dims[0]*dims[1];

  for (int i = 0; i < mask_idx.length(); i++) {
    NumericVector ret = box_nbhd(arr, dims, cds(i,0)-1, cds(i,1)-1, cds(i,2)-1, window, local, slicedim);
    NumericVector bilateral_wts = wts(_, 0) * exp(-pow(ret - arr[mask_idx[i]-1], 2) / wts(_, 1));
    double total_weight = sum(bilateral_wts);
    out[mask_idx[i]-1] = sum(ret * bilateral_wts) / total_weight;
  }

  out.attr("dim") = dims;
  return out;
}

NumericMatrix bilateral_weights_4d(int spatial_window, int temporal_window, double spatial_sigma,
                                   double intensity_sigma, double temporal_sigma, NumericVector spacing) {
  int total_elements = pow((spatial_window * 2) + 1, 3) * ((temporal_window * 2) + 1);
  NumericMatrix out(total_elements, 2);

  double spatial_denom = 2 * pow(spatial_sigma, 2);
  double intensity_denom = 2 * pow(intensity_sigma, 2);
  double temporal_denom = 2 * pow(temporal_sigma, 2);

  int ind = 0;
  for (int t = -temporal_window; t <= temporal_window; t++) {
    for (int i = -spatial_window; i <= spatial_window; i++) {
      for (int j = -spatial_window; j <= spatial_window; j++) {
        for (int k = -spatial_window; k <= spatial_window; k++) {
          out(ind, 0) = std::exp(-std::pow(i * spacing[0], 2) / spatial_denom) * std::exp(-std::pow(j * spacing[1], 2) / spatial_denom) * std::exp(-std::pow(k * spacing[2], 2) / spatial_denom) * std::exp(-std::pow(t * spacing[3], 2) / temporal_denom);
          out(ind, 1) = intensity_denom;
          ind = ind + 1;
        }
      }
    }
  }

  return out;
}


// [[Rcpp::export]]
int gridToIndexSingleCpp(IntegerVector coords, IntegerVector array_dim) {
  int rank = array_dim.size();
  int index = coords[0] - 1;

  for (int i = 1; i < rank; ++i) {
    index += (coords[i] - 1) * array_dim[i - 1];
  }

  return index;
}

// Helper function to convert 1D index to 4D index
inline IntegerVector get_4d_idx(int idx, IntegerVector& dims) {
  IntegerVector idx_4d(4);
  idx_4d[0] = idx % dims[0];
  idx_4d[1] = (idx / dims[0]) % dims[1];
  idx_4d[2] = (idx / (dims[0] * dims[1])) % dims[2];
  idx_4d[3] = idx / (dims[0] * dims[1] * dims[2]);
  return idx_4d;
}

// Helper function to convert 4D index to 1D index
inline int get_1d_idx(IntegerVector& idx_4d, IntegerVector& dims) {
  return idx_4d[0] + idx_4d[1] * dims[0] + idx_4d[2] * dims[0] * dims[1] + idx_4d[3] * dims[0] * dims[1] * dims[2];
}

inline IntegerVector get_3d_idx(int idx, IntegerVector& dims) {
  IntegerVector idx_3d(3);
  idx_3d[0] = idx % dims[0];
  idx_3d[1] = (idx / dims[0]) % dims[1];
  idx_3d[2] = idx / (dims[0] * dims[1]);
  return idx_3d;
}


// [[Rcpp::export]]
NumericVector bilateral_filter_4d_cpp(NumericVector arr, IntegerVector mask_idx, int spatial_window, int temporal_window,
                                      double spatial_sigma, double intensity_sigma, double temporal_sigma, double intensity_sd, NumericVector spacing) {
  // Extract dimensions of the input 4D array
  IntegerVector dims = arr.attr("dim");

  // Initialize the output array
  NumericVector output(arr.size());
  std::copy(arr.begin(), arr.end(), output.begin());

  // Iterate through each voxel in the mask_idx
  for (int m = 0; m < mask_idx.size(); ++m) {
    int spatial_idx = mask_idx[m];
    IntegerVector spatial_idx_3d = get_3d_idx(spatial_idx, dims);

    for (int t = 0; t < dims[3]; ++t) {
      IntegerVector idx_4d = spatial_idx_3d;
      idx_4d.push_back(t);  // Adding the time dimension to the spatial_idx_3d

      int idx = get_1d_idx(idx_4d, dims);

      double weight_sum = 0.0;
      double weighted_value_sum = 0.0;

      // Iterate through spatial and temporal neighbors
      for (int t_offset = -temporal_window; t_offset <= temporal_window; ++t_offset) {
        for (int z_offset = -spatial_window; z_offset <= spatial_window; ++z_offset) {
          for (int y_offset = -spatial_window; y_offset <= spatial_window; ++y_offset) {
            for (int x_offset = -spatial_window; x_offset <= spatial_window; ++x_offset) {
              IntegerVector neighbor_idx_4d = idx_4d + IntegerVector::create(x_offset, y_offset, z_offset, t_offset);

              // Check if the neighbor is within the bounds of the array
              if (neighbor_idx_4d[0] >= 0 && neighbor_idx_4d[0] < dims[0] &&
                  neighbor_idx_4d[1] >= 0 && neighbor_idx_4d[1] < dims[1] &&
                  neighbor_idx_4d[2] >= 0 && neighbor_idx_4d[2] < dims[2] &&
                  neighbor_idx_4d[3] >= 0 && neighbor_idx_4d[3] < dims[3]) {

                int neighbor_idx = get_1d_idx(neighbor_idx_4d, dims);
                double intensity_diff = (arr[idx] - arr[neighbor_idx]) / intensity_sd;
                double spatial_diff = std::sqrt(
                  std::pow((idx_4d[0] - neighbor_idx_4d[0]) * spacing[0], 2) +
                    std::pow((idx_4d[1] - neighbor_idx_4d[1]) * spacing[1], 2) +
                    std::pow((idx_4d[2] - neighbor_idx_4d[2]) * spacing[2], 2));
                double temporal_diff = std::abs(idx_4d[3] - neighbor_idx_4d[3]) * spacing[3];

                // Calculate the weights
                double spatial_weight = std::exp(-0.5 * std::pow(spatial_diff / spatial_sigma, 2));
                double intensity_weight = std::exp(-0.5 * std::pow(intensity_diff / intensity_sigma, 2));
                double temporal_weight = std::exp(-0.5 * std::pow(temporal_diff / temporal_sigma, 2));

                // Calculate the final weight
                double final_weight = spatial_weight * intensity_weight * temporal_weight;

                weight_sum += final_weight;
                weighted_value_sum += final_weight * arr[neighbor_idx];
              }
            }
          }
        }
      }

      // Calculate and set the bilateral filter value for the current voxel
      if (weight_sum > 0) {
        output[idx] = weighted_value_sum / weight_sum;
      }
    }
  }

  output.attr("dim") = dims;
  return output;
}






// NumericVector extract_patch(NumericVector arr, IntegerVector dims, int x, int y, int z, int patch_radius) {
//   int patch_size = 2 * patch_radius + 1;
//   NumericVector patch = NumericVector(patch_size * patch_size * patch_size);
//   int index = 0;
//
//   for (int i = -patch_radius; i <= patch_radius; i++) {
//     for (int j = -patch_radius; j <= patch_radius; j++) {
//       for (int k = -patch_radius; k <= patch_radius; k++) {
//         int xi = x + i;
//         int yi = y + j;
//         int zi = z + k;
//
//         if (xi < 0 || yi < 0 || zi < 0 || xi >= dims[0] || yi >= dims[1] || zi >= dims[2]) {
//           patch[index] = 0;
//         } else {
//           patch[index] = arr[xi + yi * dims[0] + zi * dims[0] * dims[1]];
//         }
//         index++;
//       }
//     }
//   }
//
//   return patch;
// }





// [[Rcpp::export]]
NumericVector find_seqnum(NumericVector clens, NumericVector idx) {
  NumericVector out = NumericVector(idx.length());
  long long maxid = max(idx);

  for (int i=0; i < out.length(); i++) {
    long long min_idx = 0;
    long long min_val = maxid;
    for (int j = 0; j<clens.length(); j++) {
      long long delta = idx[i] - clens[j];
      if (delta >= 0 && delta < min_val) {
        min_val = delta;
        min_idx = j + 1;
      }
    }

    out[i] = min_idx;

  }

  return out;

}

// .gridToIndex <- function(dimensions, vmat) {
// D <- Reduce("*", dimensions, accumulate=TRUE)
//   apply(vmat, 1, function(vox) {
//     sum(map_dbl(length(D):2, function(i) {
//       D[i-1]*(vox[i]-1)
//     })) + vox[1]
//   })
//
// }

// [[Rcpp::export]]
long long grid_to_intvec(IntegerVector D, IntegerVector vox) {
  long long ind = 0;
  for (int j=D.length()-1; j>0; j--) {
    ind = ind + D[j-1] * (vox[j]-1);
  }
  ind = ind + vox[0];
  return ind;
}

// [[Rcpp::export]]
NumericVector exgridToIndex3DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind,
                                 IntegerVector kind) {

  IntegerVector D = IntegerVector(array_dim.length());
  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;
  }

  int nels = iind.length()*jind.length()*kind.length();
  NumericVector out = NumericVector(nels);

  int count = 0;
  for (int k= 0; k < kind.length(); k++) {
    for (int j=0; j<jind.length(); j++) {
      for (int i=0; i<iind.length(); i++) {
        out[count] = grid_to_intvec(D, IntegerVector::create(iind[i],jind[j],kind[k]));
        count++;
      }
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericVector exgridToIndex4DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind,
                                 IntegerVector kind, IntegerVector mind) {

  IntegerVector D = IntegerVector(array_dim.length());
  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;

  }

  long long nels = iind.length()*jind.length()*kind.length()*mind.length();
  NumericVector out = NumericVector(nels);

  int count = 0;
  for (int m = 0; m<mind.length(); m++) {
    for (int k= 0; k < kind.length(); k++) {
      for (int j=0; j<jind.length(); j++) {
        for (int i=0; i<iind.length(); i++) {
          out[count] = grid_to_intvec(D, IntegerVector::create(iind[i],jind[j],kind[k],mind[m]));
          count++;
        }
      }
    }
  }

  return out;
}



// [[Rcpp::export]]
IntegerVector gridToIndexCpp(IntegerVector array_dim, IntegerMatrix voxmat) {
  IntegerVector D = IntegerVector(array_dim.length());
  IntegerVector out = IntegerVector(voxmat.nrow());

  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;

  }

  for (int i=0; i < voxmat.nrow(); i++) {
    //int ind = 0;
    //for (int j=D.length()-1; j>0; j--) {
    //  ind = ind + D[j-1] * (voxmat(i,j)-1);
    //}
    //out[i] = ind + voxmat(i,0);
    IntegerVector v = voxmat(i,_ );
    out[i] = grid_to_intvec(D, v);
  }

  return out;

}




// [[Rcpp::export]]
IntegerVector gridToIndex3DCpp(IntegerVector array_dim, NumericMatrix voxmat) {
  int slicedim = array_dim[0]*array_dim[1];
  IntegerVector out = IntegerVector(voxmat.nrow());

  for (int i=0; i < voxmat.nrow(); i++) {
    out[i] = (int)((slicedim * (voxmat(i,2) -1)) + ((voxmat(i,1)-1) * array_dim(0)) + voxmat(i,0));
  }

  return out;

}

