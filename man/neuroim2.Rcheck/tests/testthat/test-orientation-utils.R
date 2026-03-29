context("orientation utils")

library(neuroim2)

test_that("affine_to_orientation returns identity orientation", {
  aff <- diag(4)
  ornt <- affine_to_orientation(aff)
  expected <- cbind(axis = c(1, 2, 3), flip = c(1, 1, 1))
  expect_equal(ornt, expected)
})

test_that("axis-code conversions round-trip", {
  start <- c("R", "A", "S")
  ornt <- axcodes_to_orientation(start)
  back <- orientation_to_axcodes(ornt)

  expect_equal(ornt, cbind(axis = c(1, 2, 3), flip = c(1, 1, 1)))
  expect_identical(back, start)
})

test_that("orientation_transform composes and apply_orientation matches expected array operations", {
  start <- axcodes_to_orientation(c("R", "A", "S"))
  end <- axcodes_to_orientation(c("A", "R", "S"))
  tx <- orientation_transform(start, end)

  arr <- array(seq_len(2 * 3 * 4), dim = c(2, 3, 4))
  out <- apply_orientation(arr, tx)
  expected <- aperm(arr, c(2, 1, 3))

  expect_equal(tx, cbind(axis = c(2, 1, 3), flip = c(1, 1, 1)))
  expect_equal(out, expected)
})

test_that("orientation_inverse_affine handles axis flips", {
  ornt <- cbind(axis = c(1, 2, 3), flip = c(-1, 1, 1))
  aff <- orientation_inverse_affine(ornt, shape = c(5, 4, 3))

  expected <- diag(4)
  expected[1, 1] <- -1
  expected[1, 4] <- 4
  expect_equal(aff, expected)
})

test_that("findAnatomy remains compatible with identity matrix", {
  orient <- neuroim2:::findAnatomy(diag(3))
  expect_identical(orient@i, LEFT_RIGHT)
  expect_identical(orient@j, POST_ANT)
  expect_identical(orient@k, INF_SUP)
})

test_that("reorient updates NeuroSpace axes to requested orientation", {
  sp <- NeuroSpace(c(6L, 7L, 8L), spacing = c(1, 1, 1))
  sp_ras <- reorient(sp, c("R", "A", "S"))

  expect_identical(axes(sp_ras)@i, RIGHT_LEFT)
  expect_identical(axes(sp_ras)@j, ANT_POST)
  expect_identical(axes(sp_ras)@k, SUP_INF)
})
