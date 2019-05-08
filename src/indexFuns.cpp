#include <Rcpp.h>
using namespace Rcpp;




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
int grid_to_intvec(IntegerVector D, IntegerVector vox) {
  int ind = 0;
  for (int j=D.length()-1; j>0; j--) {
    ind = ind + D[j-1] * (vox[j]-1);
  }
  ind = ind + vox[0];
  return ind;
}

// [[Rcpp::export]]
IntegerVector exgridToIndex3DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind,
                                 IntegerVector kind) {

  IntegerVector D = IntegerVector(array_dim.length());
  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;
  }

  int nels = iind.length()*jind.length()*kind.length();
  IntegerVector out = IntegerVector(nels);

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
IntegerVector exgridToIndex4DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind,
                                 IntegerVector kind, IntegerVector mind) {

  IntegerVector D = IntegerVector(array_dim.length());
  int cum = 1;
  for (int i = 0; i < D.length(); i++) {
    cum = cum * array_dim[i];
    D[i] = cum;

  }

  int nels = iind.length()*jind.length()*kind.length()*mind.length();
  IntegerVector out = IntegerVector(nels);

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

