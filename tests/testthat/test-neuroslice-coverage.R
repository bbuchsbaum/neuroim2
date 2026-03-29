library(testthat)
library(neuroim2)

# ---------------------------------------------------------------------------
# NeuroSlice constructor — additional cases not in test-neuroslice.R
# ---------------------------------------------------------------------------

test_that("NeuroSlice constructor from vector (non-matrix) reshapes correctly", {
  sp <- NeuroSpace(c(5L, 5L))
  vec <- rnorm(25)
  sl <- NeuroSlice(vec, sp)
  expect_s4_class(sl, "NeuroSlice")
  expect_equal(dim(sl), c(5L, 5L))
})

test_that("NeuroSlice constructor with indices creates object with correct class and dims", {
  sp <- NeuroSpace(c(10L, 10L))
  idx <- c(1L, 50L, 100L)
  vals <- c(1.0, 2.0, 3.0)
  sl <- NeuroSlice(vals, sp, indices = idx)
  expect_s4_class(sl, "NeuroSlice")
  expect_equal(dim(sl), c(10L, 10L))
  # Indexed positions should hold the supplied values
  flat <- as.vector(sl@.Data)
  expect_equal(flat[idx], vals)
})

test_that("NeuroSlice constructor errors on 3D space", {
  sp3 <- NeuroSpace(c(5L, 5L, 5L))
  expect_error(NeuroSlice(matrix(0, 5, 5), sp3))
})

test_that("NeuroSlice constructor errors on mismatched vector length", {
  sp <- NeuroSpace(c(5L, 5L))
  expect_error(NeuroSlice(rnorm(10), sp))
})

test_that("NeuroSlice constructor errors on mismatched matrix dims", {
  sp <- NeuroSpace(c(5L, 5L))
  expect_error(NeuroSlice(matrix(0, 4, 5), sp))
})

# ---------------------------------------------------------------------------
# show() method
# ---------------------------------------------------------------------------

test_that("show method for NeuroSlice produces output", {
  sp <- NeuroSpace(c(8L, 8L), spacing = c(2, 2))
  sl <- NeuroSlice(matrix(rnorm(64), 8, 8), sp)
  out <- capture.output(show(sl))
  expect_true(length(out) > 0)
  # Should mention dimensions
  expect_true(any(grepl("8", out)))
})

test_that("show method for NeuroSlice mentions range", {
  sp <- NeuroSpace(c(5L, 5L))
  sl <- NeuroSlice(matrix(seq(1, 25), 5, 5), sp)
  out <- capture.output(show(sl))
  expect_true(any(grepl("[Rr]ange|1|25", out)))
})

test_that("show method handles NAs in data without error", {
  sp <- NeuroSpace(c(5L, 5L))
  mat <- matrix(rnorm(25), 5, 5)
  mat[1, 1] <- NA
  sl <- NeuroSlice(mat, sp)
  expect_silent(capture.output(show(sl)))
})

# ---------------------------------------------------------------------------
# mapToColors — cases not covered by existing tests
# ---------------------------------------------------------------------------

test_that("mapToColors with alpha < 1 returns 4-column RGBA array", {
  vals <- matrix(runif(16), 4, 4)
  result <- mapToColors(vals, col = heat.colors(64), alpha = 0.5,
                        irange = c(0, 1))
  # Should return an array with last dim = 4 (RGBA)
  expect_equal(length(dim(result)), 3)
  expect_equal(dim(result)[3], 4)
})

test_that("mapToColors with vector input and alpha < 1 returns Nx4 matrix", {
  vals <- runif(10)
  result <- mapToColors(vals, col = heat.colors(64), alpha = 0.5,
                        irange = c(0, 1))
  expect_equal(dim(result), c(10L, 4L))
})

test_that("mapToColors errors when irange is decreasing", {
  expect_error(mapToColors(1:5, irange = c(5, 1)))
})

test_that("mapToColors clips values below irange min", {
  vals <- c(-10, 0.5, 1)
  col <- heat.colors(128)
  result <- mapToColors(vals, col = col, irange = c(0, 1), alpha = 1)
  # Value -10 should be clipped to the first color (index 1)
  expect_equal(result[1], col[1])
})

test_that("mapToColors clips values above irange max", {
  vals <- c(0.5, 1, 10)
  col <- heat.colors(128)
  result <- mapToColors(vals, col = col, irange = c(0, 1), alpha = 1)
  # Value 10 should be clipped to the last color
  expect_equal(result[3], col[length(col)])
})

# ---------------------------------------------------------------------------
# mask() method on NeuroSlice
# ---------------------------------------------------------------------------

test_that("mask() on NeuroSlice returns LogicalNeuroVol", {
  sp <- NeuroSpace(c(6L, 6L))
  sl <- NeuroSlice(matrix(rnorm(36), 6, 6), sp)
  m <- mask(sl)
  expect_s4_class(m, "LogicalNeuroVol")
})

test_that("mask() on NeuroSlice has 3 dimensions", {
  sp <- NeuroSpace(c(6L, 6L))
  sl <- NeuroSlice(matrix(rnorm(36), 6, 6), sp)
  m <- mask(sl)
  expect_equal(length(dim(m)), 3)
  # Third dim should be 1 (expanded from 2D)
  expect_equal(dim(m)[3], 1L)
})

# ---------------------------------------------------------------------------
# grid_to_index and index_to_grid — edge cases
# ---------------------------------------------------------------------------

test_that("grid_to_index on NeuroSlice with single numeric vector works", {
  sp <- NeuroSpace(c(10L, 10L))
  sl <- NeuroSlice(matrix(seq_len(100), 10, 10), sp)
  idx <- grid_to_index(sl, c(3L, 4L))
  # column-major: (row-1) + (col-1)*nrow + 1
  expect_equal(idx, as.integer((4 - 1) * 10 + 3))
})

test_that("index_to_grid round-trips with grid_to_index", {
  sp <- NeuroSpace(c(8L, 8L))
  sl <- NeuroSlice(matrix(seq_len(64), 8, 8), sp)
  coords_in <- matrix(c(2L, 3L, 5L, 7L), nrow = 2, byrow = TRUE)
  idx <- grid_to_index(sl, coords_in)
  coords_out <- index_to_grid(sl, idx)
  expect_equal(coords_out, coords_in)
})
