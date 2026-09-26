test_that("explicit sparse orientation preserves non-symmetric square series", {
  mask <- array(seq_len(8) %in% c(1, 3, 6), c(2, 2, 2))
  ix <- c(1L, 3L, 6L)

  # Rectangular: auto inference remains unambiguous either way.
  sp_rect <- NeuroSpace(c(2, 2, 2, 4))
  y_rect <- matrix(seq_len(4 * 3), 4, 3)
  expect_equal(
    series(SparseNeuroVec(y_rect, sp_rect, mask), ix),
    y_rect
  )
  expect_equal(
    series(SparseNeuroVec(t(y_rect), sp_rect, mask), ix),
    y_rect
  )

  # Square: time-by-voxels needs an explicit orientation (issue #31).
  sp <- NeuroSpace(c(2, 2, 2, 3))
  y <- matrix(seq_len(3 * 3), 3, 3)
  a <- SparseNeuroVec(y, sp, mask, orientation = "time_x_voxels")
  b <- SparseNeuroVec(t(y), sp, mask, orientation = "voxels_x_time")
  expect_equal(series(a, ix), y)
  expect_equal(series(b, ix), y)

  full <- matrix(0, 8, 3)
  full[ix, ] <- t(y)
  expect_equal(as.matrix(a), full)
  expect_equal(as.matrix(b), full)

  # Auto keeps the historic voxels-by-time convention and warns when square.
  expect_warning(
    auto_vox <- SparseNeuroVec(t(y), sp, mask),
    "assuming.*voxels_x_time"
  )
  expect_equal(series(auto_vox, ix), y)

  expect_warning(
    auto_time <- SparseNeuroVec(y, sp, mask),
    "assuming.*voxels_x_time"
  )
  # Without an explicit orientation, a square series() matrix is transposed.
  expect_false(isTRUE(all.equal(series(auto_time, ix), y)))
})

test_that("explicit sparse orientation rejects mismatched shapes and bad values", {
  mask <- array(seq_len(8) %in% c(1, 3, 6), c(2, 2, 2))
  sp <- NeuroSpace(c(2, 2, 2, 4))

  expect_error(
    SparseNeuroVec(matrix(1, 4, 3), sp, mask, orientation = "voxels_x_time"),
    "orientation"
  )
  expect_error(
    SparseNeuroVec(matrix(1, 4, 3), sp, mask, orientation = "bad"),
    "arg"
  )
})
