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
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

struct RowRepresentativeWorker : public Worker {
  const RMatrix<double> mat;
  RVector<double> out;
  const bool use_median;

  RowRepresentativeWorker(const NumericMatrix& mat,
                          NumericVector& out,
                          bool use_median)
    : mat(mat), out(out), use_median(use_median) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> buffer;
    if (use_median) {
      buffer.reserve(mat.ncol());
    }

    for (std::size_t ii = begin; ii < end; ++ii) {
      if (use_median) {
        buffer.clear();
        bool has_missing = false;

        for (int jj = 0; jj < mat.ncol(); ++jj) {
          double value = mat(ii, jj);
          if (NumericVector::is_na(value) || ISNAN(value)) {
            has_missing = true;
            break;
          }
          buffer.push_back(value);
        }

        if (has_missing || buffer.empty()) {
          out[ii] = NA_REAL;
          continue;
        }

        std::size_t mid = buffer.size() / 2;
        std::nth_element(buffer.begin(), buffer.begin() + mid, buffer.end());
        double upper = buffer[mid];

        if ((buffer.size() % 2U) == 1U) {
          out[ii] = upper;
        } else {
          std::nth_element(buffer.begin(), buffer.begin() + mid - 1, buffer.begin() + mid);
          double lower = buffer[mid - 1];
          out[ii] = (lower + upper) / 2.0;
        }
      } else {
        long double sum = 0.0L;
        bool has_missing = false;

        for (int jj = 0; jj < mat.ncol(); ++jj) {
          double value = mat(ii, jj);
          if (NumericVector::is_na(value) || ISNAN(value)) {
            has_missing = true;
            break;
          }
          sum += std::fabs(value);
        }

        out[ii] = has_missing ? NA_REAL : static_cast<double>(sum / static_cast<long double>(mat.ncol()));
      }
    }
  }
};

}  // namespace

// [[Rcpp::export]]
NumericVector representative_volume_cpp(NumericMatrix mat, std::string representative) {
  if (mat.nrow() == 0) {
    return NumericVector(0);
  }

  if (mat.ncol() == 0) {
    stop("mat must have at least one column");
  }

  bool use_median;
  if (representative == "median") {
    use_median = true;
  } else if (representative == "mean_abs") {
    use_median = false;
  } else {
    stop("representative must be one of 'median' or 'mean_abs'");
  }

  NumericVector out(mat.nrow());
  RowRepresentativeWorker worker(mat, out, use_median);
  parallelFor(0, static_cast<std::size_t>(mat.nrow()), worker);
  return out;
}
