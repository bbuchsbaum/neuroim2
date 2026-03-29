context("summary methods for NeuroVol and NeuroVec")

# ---- DenseNeuroVol ----

test_that("summary(DenseNeuroVol) returns class summary.NeuroVol", {
  vol <- make_vol(dim = c(4L, 4L, 4L))
  s <- summary(vol)
  expect_s3_class(s, "summary.NeuroVol")
})

test_that("summary(DenseNeuroVol) has correct required fields", {
  vol <- make_vol(dim = c(4L, 4L, 4L))
  s <- summary(vol)
  expected_fields <- c("class", "dim", "spacing", "origin", "orientation",
                       "min", "max", "mean", "sd", "NAs", "zeros", "nonzeros")
  expect_true(all(expected_fields %in% names(s)))
})

test_that("summary(DenseNeuroVol) dim matches object dim", {
  vol <- make_vol(dim = c(5L, 6L, 7L))
  s <- summary(vol)
  expect_equal(s$dim, c(5L, 6L, 7L))
})

test_that("summary(DenseNeuroVol) class field is correct string", {
  vol <- make_vol()
  s <- summary(vol)
  expect_equal(s$class, "DenseNeuroVol")
})

test_that("summary(DenseNeuroVol) min/max/mean are numeric scalars", {
  set.seed(1)
  vol <- make_vol(dim = c(4L, 4L, 4L), data = 1:64)
  s <- summary(vol)
  expect_equal(s$min, 1)
  expect_equal(s$max, 64)
  expect_equal(s$mean, mean(1:64))
})

test_that("summary(DenseNeuroVol) NAs field counts NAs correctly", {
  d <- array(1:64, c(4L, 4L, 4L))
  d[1, 1, 1] <- NA
  d[2, 2, 2] <- NA
  vol <- DenseNeuroVol(d, NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1)))
  s <- summary(vol)
  expect_equal(s$NAs, 2L)
})

test_that("summary(DenseNeuroVol) zeros field counts zeros correctly", {
  d <- array(c(rep(0, 10), rep(1, 54)), c(4L, 4L, 4L))
  vol <- DenseNeuroVol(d, NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1)))
  s <- summary(vol)
  expect_equal(s$zeros, 10L)
  expect_equal(s$nonzeros, 54L)
})

test_that("summary(DenseNeuroVol) all-zero volume", {
  vol <- DenseNeuroVol(array(0, c(3L, 3L, 3L)),
                       NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1)))
  s <- summary(vol)
  expect_equal(s$min, 0)
  expect_equal(s$max, 0)
  expect_equal(s$zeros, 27L)
  expect_equal(s$nonzeros, 0L)
})

test_that("print.summary.NeuroVol produces output", {
  vol <- make_vol(dim = c(4L, 4L, 4L))
  s <- summary(vol)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("DenseNeuroVol", out)))
})

test_that("print.summary.NeuroVol shows NA line only when NAs > 0", {
  vol_clean <- make_vol(dim = c(3L, 3L, 3L), data = rep(1, 27))
  s_clean <- summary(vol_clean)
  out_clean <- capture.output(print(s_clean))
  expect_false(any(grepl("NAs", out_clean)))

  d <- array(1, c(3L, 3L, 3L))
  d[1] <- NA
  vol_na <- DenseNeuroVol(d, NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1)))
  s_na <- summary(vol_na)
  out_na <- capture.output(print(s_na))
  expect_true(any(grepl("NAs", out_na)))
})

# ---- DenseNeuroVec ----

test_that("summary(DenseNeuroVec) returns class summary.NeuroVec", {
  vec <- make_vec(dim = c(4L, 4L, 4L), ntime = 10L)
  s <- summary(vec)
  expect_s3_class(s, "summary.NeuroVec")
})

test_that("summary(DenseNeuroVec) has correct required fields", {
  vec <- make_vec(dim = c(4L, 4L, 4L), ntime = 10L)
  s <- summary(vec)
  expected_fields <- c("class", "dim", "spacing", "origin", "orientation",
                       "time_points", "global_min", "global_max", "global_mean",
                       "temporal_mean_range", "temporal_sd_range",
                       "NAs", "nonzero_voxels", "total_voxels")
  expect_true(all(expected_fields %in% names(s)))
})

test_that("summary(DenseNeuroVec) time_points matches 4th dim", {
  vec <- make_vec(dim = c(4L, 4L, 4L), ntime = 15L)
  s <- summary(vec)
  expect_equal(s$time_points, 15L)
})

test_that("summary(DenseNeuroVec) total_voxels matches spatial dims", {
  vec <- make_vec(dim = c(4L, 5L, 6L), ntime = 8L)
  s <- summary(vec)
  expect_equal(s$total_voxels, 4L * 5L * 6L)
})

test_that("summary(DenseNeuroVec) global_mean matches manual calculation", {
  set.seed(42)
  dat <- array(rnorm(4 * 4 * 4 * 10), c(4L, 4L, 4L, 10L))
  sp <- NeuroSpace(c(4L, 4L, 4L, 10L), c(1, 1, 1))
  vec <- DenseNeuroVec(dat, sp)
  s <- summary(vec)
  expect_equal(s$global_mean, mean(dat), tolerance = 1e-10)
})

test_that("print.summary.NeuroVec produces output", {
  vec <- make_vec(dim = c(4L, 4L, 4L), ntime = 10L)
  s <- summary(vec)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("DenseNeuroVec", out)))
  expect_true(any(grepl("Time points", out)))
})

# ---- SparseNeuroVec ----

test_that("summary(SparseNeuroVec) returns class summary.NeuroVec", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  s <- summary(svec)
  expect_s3_class(s, "summary.NeuroVec")
})

test_that("summary(SparseNeuroVec) nonzero_voxels matches indices length", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  s <- summary(svec)
  expect_equal(s$nonzero_voxels, length(indices(svec)))
})

test_that("summary(SparseNeuroVec) total_voxels matches spatial dims", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  s <- summary(svec)
  expect_equal(s$total_voxels, 5L * 5L * 5L)
})

test_that("summary(SparseNeuroVec) time_points correct", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 12L, frac = 0.3)
  s <- summary(svec)
  expect_equal(s$time_points, 12L)
})

test_that("print.summary.NeuroVec works for SparseNeuroVec summary", {
  svec <- make_sparse_vec(dim = c(5L, 5L, 5L), ntime = 8L, frac = 0.4)
  s <- summary(svec)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("SparseNeuroVec", out)))
})
