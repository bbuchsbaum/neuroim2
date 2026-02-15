context("affine utils")

library(neuroim2)

test_that("apply_affine transforms matrix and preserves array leading dims", {
  aff <- diag(c(2, 3, 4, 1))
  aff[1:3, 4] <- c(10, 20, 30)

  pts <- rbind(c(1, 2, 3), c(4, 5, 6))
  out <- apply_affine(aff, pts)
  expected <- pts %*% t(aff[1:3, 1:3]) + matrix(aff[1:3, 4], nrow = nrow(pts), ncol = 3, byrow = TRUE)
  expect_equal(out, expected)

  arr <- array(pts, dim = c(1, 2, 3))
  out_arr <- apply_affine(aff, arr)
  expect_equal(dim(out_arr), c(1, 2, 3))
  expect_equal(as.numeric(out_arr), as.numeric(expected))
})

test_that("to_matvec and from_matvec roundtrip", {
  aff <- diag(c(2, 3, 4, 1))
  aff[1:3, 4] <- c(9, 10, 11)

  mv <- to_matvec(aff)
  aff2 <- from_matvec(mv$matrix, mv$vector)

  expect_equal(mv$matrix, aff[1:3, 1:3])
  expect_equal(mv$vector, aff[1:3, 4])
  expect_equal(aff2, aff)
})

test_that("append_diag expands affine and validates starts", {
  aff <- diag(4)
  out <- append_diag(aff, steps = c(9, 10), starts = c(99, 100))

  expect_equal(dim(out), c(6, 6))
  expect_equal(diag(out)[4:5], c(9, 10))
  expect_equal(out[4:5, 6], c(99, 100))
  expect_equal(out[6, 6], 1)

  expect_error(append_diag(aff, steps = c(1, 2), starts = 1), class = "AffineError")
})

test_that("dot_reduce computes right-associated matrix product", {
  A <- matrix(c(1, 2, 3, 4), 2, 2)
  B <- matrix(c(0, 1, 1, 0), 2, 2)
  C <- matrix(c(2, 0, 0, 2), 2, 2)

  expect_equal(dot_reduce(A, B, C), A %*% (B %*% C))
})

test_that("voxel_sizes and obliquity are consistent for cardinal axes", {
  aff <- diag(c(2, 3, 4, 1))
  expect_equal(voxel_sizes(aff), c(2, 3, 4))
  expect_equal(obliquity(aff), c(0, 0, 0))
})

test_that("rescale_affine preserves center world coordinate", {
  aff <- diag(c(2, 3, 4, 1))
  aff[1:3, 4] <- c(10, 20, 30)

  shape <- c(5, 7, 9)
  new_shape <- c(9, 11, 13)
  new_zooms <- c(1, 1.5, 2)

  out <- rescale_affine(aff, shape = shape, zooms = new_zooms, new_shape = new_shape)
  expect_equal(voxel_sizes(out), new_zooms)

  c_in <- floor((shape - 1) / 2)
  c_out <- floor((new_shape - 1) / 2)
  expect_equal(apply_affine(aff, c_in), apply_affine(out, c_out))
})
