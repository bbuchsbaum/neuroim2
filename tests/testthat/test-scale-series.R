test_that("scale_series(SparseNeuroVec) matches the historical column-wise formula", {
  set.seed(1)
  arr <- array(rnorm(6 * 5 * 4 * 7), dim = c(6, 5, 4, 7))
  dense <- DenseNeuroVec(arr, NeuroSpace(dim(arr)))
  mask_arr <- array(runif(prod(dim(arr)[1:3])) > 0.65, dim(arr)[1:3])
  sparse <- as.sparse(dense, LogicalNeuroVol(mask_arr, drop_dim(space(dense))))

  ref_scale <- function(M, center, scale) {
    if (center) {
      M <- sweep(M, 2, colMeans(M), "-")
    }
    if (scale) {
      if (nrow(M) <= 1L) {
        sds <- rep(1, ncol(M))
      } else if (center) {
        sds <- sqrt(colSums(M * M) / (nrow(M) - 1L))
      } else {
        centered <- sweep(M, 2, colMeans(M), "-")
        sds <- sqrt(colSums(centered * centered) / (nrow(M) - 1L))
      }
      sds[!is.finite(sds) | sds == 0] <- 1
      M <- sweep(M, 2, sds, "/")
    }
    unname(as.matrix(M))
  }

  scaled_only <- scale_series(sparse, center = FALSE, scale = TRUE)
  centered_only <- scale_series(sparse, center = TRUE, scale = FALSE)
  centered_scaled <- scale_series(sparse, center = TRUE, scale = TRUE)

  expect_equal(scaled_only@data, ref_scale(sparse@data, FALSE, TRUE))
  expect_equal(centered_only@data, ref_scale(sparse@data, TRUE, FALSE))
  expect_equal(centered_scaled@data, ref_scale(sparse@data, TRUE, TRUE))
})
