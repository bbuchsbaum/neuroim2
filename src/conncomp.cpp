#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif
using namespace Rcpp;



// neighbors <- function(vox) {
//
//   vox.hood <- t(tlocal.mask + vox)
//   if (any(vox == 1) || any(vox == DIM)) {
//     vox.hood <- vox.hood[apply(vox.hood, 1, function(coords) {
//       all(coords > 1 & coords <= DIM)
//     }),,drop=FALSE]
//   }
//
//   vox.hood[labels[vox.hood] != 0,,drop=F]
// }





// int find(IntegerVector nodes, int i) {
//
//   while (nodes(i) != i) {
//     i = nodes(i);
//   }
//
//   return nodes(i);
// }
//
//
// NumericVector neighbors(NumericVector vox, NumericMatrix local_mask, IntegerVector dim) {
//   NumericMatrix vox_hood(local_mask.nrow(), local_mask.ncol());
//   LogicalVector keep(local_mask.nrow(), true);
//
//   int count = 0;
//   for (int i=0; i<local_mask.nrow(); i++) {
//     vox_hood(i,0) = vox[0] + local_mask(i,0);
//     vox_hood(i,1) = vox[1] + local_mask(i,1);
//     vox_hood(i,2) = vox[2] + local_mask(i,2);
//
//     bool border = false;
//     for (int j=0; j<2; j++) {
//       if ( (vox_hood(i,j) < 1) || (vox_hood(i,j) > dim(j))) {
//         border=true;
//         break;
//       }
//     }
//
//     if (border) {
//       keep(i) = false;
//     } else {
//       count++;
//     }
//   }
//
//   if (count > 0) {
//     NumericMatrix hood(count, local_mask.ncol());
//     int index=0;
//     for (int i=0; i<keep.length(); i++) {
//       if (keep(i)) {
//         hood(index, _) = vox_hood(i,_);
//         index++;
//       }
//     }
//     return hood;
//   } else {
//     return vox_hood;
//   }
// }

// IntegerVector conncomp3dCpp(NumericMatrix grid, NumericMatrix mask, IntegerVector dim) {
//   for (int i =0; i<grid.nrow(); i++) {
//     IntegerVector vox = grid(i,_);
//     NumericMatrix nabes = neighbors(vox, mask, dim);
//     if (nrow(nabes) == 0) {
//     nodes[nextlabel] <- nextlabel
//     labels[vox[1],vox[2],vox[3]] <- nextlabel
//   } else {
//     L <- labels[nabes]
//     ML <- min(L)
//     labels[vox[1],vox[2], vox[3]] <- ML
//     nodes[nextlabel] <- ML
//     for (lab in L) {
//       rootx <- find(lab)
//       nodes[rootx] <- find(ML)
//     }
//   }
//
//   nextlabel <- nextlabel + 1
// }

// pass2
// for (k in 1:zdim) {
//   for (j in 1:ydim) {
//     for (i in 1:xdim) {
//       if (labels[i,j,k] > 0) {
//         labels[i,j,k] <- find(labels[i,j,k])
//       }
//     }
//   }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
local_mask = as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1), z=c(-1,0,1)))
vox = c(10,10,10)
neighbors(vox, local_mask, dim=c(15,10,15))
*/
