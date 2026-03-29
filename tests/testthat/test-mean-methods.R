context("mean methods for NeuroVec types")

# ---- mean(DenseNeuroVec) ----

test_that("mean(DenseNeuroVec) returns a DenseNeuroVol", {
  vec <- make_vec(dim = c(5L, 5L, 5L), ntime = 10L)
  result <- mean(vec)
  expect_s4_class(result, "DenseNeuroVol")
})

test_that("mean(DenseNeuroVec) result has correct 3D dims", {
  vec <- make_vec(dim = c(4L, 5L, 6L), ntime = 8L)
  result <- mean(vec)
  expect_equal(dim(result), c(4L, 5L, 6L))
})

test_that("mean(DenseNeuroVec) of constant vec equals that constant", {
  sp4 <- NeuroSpace(c(3L, 3L, 3L, 5L), c(1, 1, 1))
  dat <- array(7, c(3L, 3L, 3L, 5L))
  vec <- DenseNeuroVec(dat, sp4)
  result <- mean(vec)
  expect_equal(as.numeric(result@.Data), rep(7, 27))
})

test_that("mean(DenseNeuroVec) matches rowMeans of matrix reshape", {
  set.seed(99)
  dim3 <- c(4L, 4L, 4L)
  ntime <- 12L
  sp4 <- NeuroSpace(c(dim3, ntime), c(1, 1, 1))
  dat <- array(rnorm(prod(dim3) * ntime), c(dim3, ntime))
  vec <- DenseNeuroVec(dat, sp4)

  result <- mean(vec)
  M <- matrix(dat, nrow = prod(dim3), ncol = ntime)
  expected <- rowMeans(M)

  expect_equal(as.numeric(result@.Data), expected, tolerance = 1e-10)
})

test_that("mean(DenseNeuroVec) space matches drop_dim of input space", {
  vec <- make_vec(dim = c(4L, 4L, 4L), ntime = 8L)
  result <- mean(vec)
  expect_equal(dim(space(result)), c(4L, 4L, 4L))
  expect_equal(spacing(space(result)), spacing(space(vec))[1:3])
})

test_that("mean(DenseNeuroVec) single time point returns zero-variance result", {
  sp4 <- NeuroSpace(c(3L, 3L, 3L, 1L), c(1, 1, 1))
  dat <- array(rnorm(27), c(3L, 3L, 3L, 1L))
  vec <- DenseNeuroVec(dat, sp4)
  result <- mean(vec)
  expect_equal(as.numeric(result@.Data), as.numeric(dat), tolerance = 1e-10)
})

# ---- mean(SparseNeuroVec) ----

test_that("mean(SparseNeuroVec) returns a SparseNeuroVol", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  result <- mean(svec)
  expect_s4_class(result, "SparseNeuroVol")
})

test_that("mean(SparseNeuroVec) result has correct 3D dims", {
  svec <- make_sparse_vec(dim = c(4L, 5L, 6L), ntime = 8L, frac = 0.4)
  result <- mean(svec)
  expect_equal(dim(result), c(4L, 5L, 6L))
})

test_that("mean(SparseNeuroVec) indices match original SparseNeuroVec", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  result <- mean(svec)
  # SparseNeuroVol has no indices() method; use @data@i for nonzero positions
  expect_equal(sort(result@data@i), sort(indices(svec)))
})

test_that("mean(SparseNeuroVec) values match colMeans of data matrix", {
  set.seed(7)
  sdim <- c(5L, 5L, 5L)
  ntime <- 10L
  sp4 <- NeuroSpace(c(sdim, ntime), c(1, 1, 1))
  mask <- array(runif(prod(sdim)) < 0.4, sdim)
  nk <- sum(mask)
  dat <- matrix(rnorm(ntime * nk), nrow = ntime, ncol = nk)
  svec <- SparseNeuroVec(dat, sp4, mask = mask)

  result <- mean(svec)
  expected_vals <- colMeans(svec@data)

  # Extract result values at the active indices in the same order
  result_vals <- as.numeric(result@data)[result@data@i]
  expect_equal(result_vals, expected_vals, tolerance = 1e-10)
})

test_that("mean(SparseNeuroVec) of constant data equals that constant", {
  sdim <- c(4L, 4L, 4L)
  ntime <- 6L
  sp4 <- NeuroSpace(c(sdim, ntime), c(1, 1, 1))
  mask <- array(TRUE, sdim)
  nk <- sum(mask)
  dat <- matrix(3, nrow = ntime, ncol = nk)
  svec <- SparseNeuroVec(dat, sp4, mask = mask)
  result <- mean(svec)
  active_vals <- as.numeric(result@data)[result@data@i]
  expect_equal(unique(active_vals), 3, tolerance = 1e-10)
})
