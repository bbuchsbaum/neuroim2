context("performance guardrails")

library(neuroim2)

skip_if_perf_disabled <- function() {
  if (!identical(Sys.getenv("NEUROIM2_PERF"), "true")) {
    skip("Performance guardrails run only when NEUROIM2_PERF=true")
  }
}

old_sparse_subset_reference <- function(x, i, j, k, m, drop = TRUE) {
  grid_to_index3d <- get(".gridToIndex3D", envir = asNamespace("neuroim2"))
  vmat <- as.matrix(expand.grid(i, j, k))
  ind <- grid_to_index3d(dim(x)[1:3], vmat[, 1:3, drop = FALSE])
  mapped <- lookup(x, ind)
  keep <- mapped > 0
  dimout <- c(length(i), length(j), length(k), length(m))

  if (!any(keep)) {
    out0 <- array(0, dimout)
    return(if (drop) base::drop(out0) else out0)
  }

  out <- array(0, dimout)
  nz_idx <- which(keep)
  vals <- matricized_access(x, mapped[nz_idx])
  vals <- vals[m, , drop = FALSE]
  coords <- vmat[nz_idx, , drop = FALSE]
  ii_pos <- match(coords[, 1], i)
  jj_pos <- match(coords[, 2], j)
  kk_pos <- match(coords[, 3], k)
  for (col in seq_along(nz_idx)) {
    out[ii_pos[col], jj_pos[col], kk_pos[col], ] <- vals[, col]
  }

  if (drop) base::drop(out) else out
}

test_that("dense matrix series fast path beats scalar extraction", {
  skip_on_cran()
  skip_if_perf_disabled()

  set.seed(1)
  dims <- c(18, 18, 18, 40)
  x <- DenseNeuroVec(array(rnorm(prod(dims)), dim = dims), NeuroSpace(dims))
  coords <- cbind(
    sample.int(dims[1], 400, replace = TRUE),
    sample.int(dims[2], 400, replace = TRUE),
    sample.int(dims[3], 400, replace = TRUE)
  )

  fast_time <- system.time(replicate(
    10,
    series(x, coords),
    simplify = FALSE
  ))[["elapsed"]]

  ref_time <- system.time(replicate(
    10,
    vapply(
      seq_len(nrow(coords)),
      function(ii) series(x, coords[ii, 1], coords[ii, 2], coords[ii, 3]),
      numeric(dims[4])
    ),
    simplify = FALSE
  ))[["elapsed"]]

  expect_lt(fast_time, ref_time)
})

test_that("sparse subset path beats the old expand.grid implementation", {
  skip_on_cran()
  skip_if_perf_disabled()

  set.seed(1)
  dims <- c(24, 24, 16, 10)
  arr <- array(rnorm(prod(dims)), dim = dims)
  dense <- DenseNeuroVec(arr, NeuroSpace(dims))
  mask_arr <- array(runif(prod(dims[1:3])) > 0.92, dims[1:3])
  sparse <- as.sparse(dense, LogicalNeuroVol(mask_arr, drop_dim(space(dense))))
  i <- 3:18
  j <- 4:17
  k <- 2:11
  m <- 1:6

  fast_time <- system.time(replicate(
    20,
    sparse[i, j, k, m, drop = FALSE],
    simplify = FALSE
  ))[["elapsed"]]

  ref_time <- system.time(replicate(
    20,
    old_sparse_subset_reference(sparse, i, j, k, m, drop = FALSE),
    simplify = FALSE
  ))[["elapsed"]]

  expect_lt(fast_time, ref_time)
})
