context("plot helpers coverage")

library(neuroim2)

# ---- slice_to_matrix ----

test_that("slice_to_matrix handles a plain matrix", {
  m <- matrix(1:9, nrow = 3)
  result <- neuroim2:::slice_to_matrix(m)
  expect_true(is.matrix(result))
  expect_equal(result, m)
})

test_that("slice_to_matrix handles a NeuroSlice", {
  sp <- NeuroSpace(c(5L, 5L), c(1, 1))
  sl <- NeuroSlice(array(1:25, c(5, 5)), sp)
  result <- neuroim2:::slice_to_matrix(sl)
  expect_true(is.matrix(result))
})

# ---- slice_df ----

test_that("slice_df returns a data.frame with x, y, value columns", {
  m <- matrix(rnorm(25), nrow = 5)
  df <- neuroim2:::slice_df(m)
  expect_true(is.data.frame(df))
  expect_true(all(c("x", "y", "value") %in% names(df)))
  expect_equal(nrow(df), 25L)
})

test_that("slice_df respects downsample parameter", {
  m <- matrix(rnorm(100), nrow = 10)
  df_full  <- neuroim2:::slice_df(m, downsample = 1L)
  df_down  <- neuroim2:::slice_df(m, downsample = 2L)
  expect_lt(nrow(df_down), nrow(df_full))
})

test_that("slice_df works with a NeuroSlice", {
  sp <- NeuroSpace(c(6L, 6L), c(1, 1))
  sl <- NeuroSlice(array(rnorm(36), c(6, 6)), sp)
  df <- neuroim2:::slice_df(sl)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 36L)
})

# ---- compute_limits ----

test_that("compute_limits mode='data' returns min/max", {
  x <- c(1, 2, 3, 4, 5)
  lim <- neuroim2:::compute_limits(x, mode = "data")
  expect_equal(lim, c(1, 5))
})

test_that("compute_limits mode='robust' returns quantile range", {
  set.seed(1)
  x <- rnorm(1000)
  lim <- neuroim2:::compute_limits(x, mode = "robust", probs = c(0.02, 0.98))
  expect_length(lim, 2)
  expect_lt(lim[1], lim[2])
})

test_that("compute_limits returns c(0,1) for empty input", {
  lim <- neuroim2:::compute_limits(numeric(0))
  expect_equal(lim, c(0, 1))
})

test_that("compute_limits handles all-NA input", {
  lim <- neuroim2:::compute_limits(c(NA_real_, NA_real_))
  expect_equal(lim, c(0, 1))
})

test_that("compute_limits falls back to data range when robust range collapses", {
  x <- rep(3, 100)
  lim <- neuroim2:::compute_limits(x, mode = "robust")
  expect_equal(lim[1], lim[2])  # range is trivially equal when all same
})

# ---- rescale01 ----

test_that("rescale01 maps to [0,1]", {
  x <- c(0, 5, 10)
  r <- neuroim2:::rescale01(x, from = c(0, 10))
  expect_equal(r, c(0, 0.5, 1))
})

test_that("rescale01 returns 0.5 when from range is zero", {
  x <- c(1, 2, 3)
  r <- neuroim2:::rescale01(x, from = c(5, 5))
  expect_true(all(r == 0.5))
})

test_that("rescale01 returns x unchanged when from is NULL", {
  x <- c(1, 2, 3)
  r <- neuroim2:::rescale01(x, from = NULL)
  expect_equal(r, x)
})

test_that("rescale01 returns x unchanged when from has non-finite values", {
  x <- c(1, 2, 3)
  r <- neuroim2:::rescale01(x, from = c(NA_real_, 5))
  expect_equal(r, x)
})

# ---- matrix_to_colors ----

test_that("matrix_to_colors returns a character vector", {
  m <- matrix(1:9, nrow = 3)
  cols <- neuroim2:::matrix_to_colors(m)
  expect_type(cols, "character")
  expect_equal(length(cols), 9L)
})

test_that("matrix_to_colors respects explicit limits", {
  m <- matrix(c(0, 5, 10), nrow = 1)
  cols <- neuroim2:::matrix_to_colors(m, limits = c(0, 10))
  expect_type(cols, "character")
  expect_equal(length(cols), 3L)
})

test_that("matrix_to_colors works with alpha_map", {
  m <- matrix(1:4, nrow = 2)
  amap <- matrix(c(0.5, 0.5, 1, 1), nrow = 2)
  cols <- neuroim2:::matrix_to_colors(m, alpha_map = amap)
  expect_type(cols, "character")
})

# ---- matrix_to_rgba ----

test_that("matrix_to_rgba returns a 4-layer array", {
  m <- matrix(1:9, nrow = 3)
  rgba <- neuroim2:::matrix_to_rgba(m)
  expect_equal(length(dim(rgba)), 3L)
  expect_equal(dim(rgba)[3], 4L)
})

test_that("matrix_to_rgba alpha channel is between 0 and 1", {
  m <- matrix(1:9, nrow = 3)
  rgba <- neuroim2:::matrix_to_rgba(m)
  expect_true(all(rgba[,,4] >= 0 & rgba[,,4] <= 1))
})

test_that("matrix_to_rgba with alpha_map overrides alpha", {
  m <- matrix(1:4, nrow = 2)
  amap <- matrix(c(0.2, 0.4, 0.6, 0.8), nrow = 2)
  rgba <- neuroim2:::matrix_to_rgba(m, alpha_map = amap)
  expect_equal(dim(rgba), c(2L, 2L, 4L))
})

test_that("matrix_to_rgba with NULL limits auto-computes", {
  m <- matrix(seq(0, 1, length.out = 9), nrow = 3)
  rgba <- neuroim2:::matrix_to_rgba(m, limits = NULL)
  expect_equal(dim(rgba)[3], 4L)
})

# ---- annotate_orientation ----

test_that("annotate_orientation returns a list of ggplot layers", {
  layers <- annotate_orientation("axial", dims = c(10, 10))
  expect_true(is.list(layers))
  expect_equal(length(layers), 4L)
})

test_that("annotate_orientation works for coronal plane", {
  layers <- annotate_orientation("coronal", dims = c(8, 12))
  expect_equal(length(layers), 4L)
})

test_that("annotate_orientation works for sagittal plane", {
  layers <- annotate_orientation("sagittal", dims = c(8, 12))
  expect_equal(length(layers), 4L)
})

# ---- coord_neuro_fixed ----

test_that("coord_neuro_fixed returns a list of ggplot components", {
  result <- neuroim2:::coord_neuro_fixed()
  expect_true(is.list(result))
  expect_equal(length(result), 2L)
})
