#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector find_seqnum(IntegerVector clens, IntegerVector idx) {
  IntegerVector out = IntegerVector(idx.length());
  int maxid = max(idx);

  for (int i=0; i < out.length(); i++) {
    int min_idx = 0;
    int min_val = maxid;
    for (int j = 0; j<clens.length(); j++) {
      int delta = idx[i] - clens[j];
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
IntegerVector gridToIndexCpp(IntegerVector array_dim, NumericMatrix voxmat) {
  IntegerVector D = IntegerVector(array_dim.length());
  IntegerVector out = IntegerVector(voxmat.nrow());

  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;

  }

  for (int i=0; i < voxmat.nrow(); i++) {
    int ind = 0;
    for (int j=D.length()-1; j>0; j--) {
      ind = ind + D[j-1] * (voxmat(i,j)-1);
    }
    out[i] = ind + voxmat(i,0);
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

