context("space utils")

library(neuroim2)

test_that("output_aligned_space returns identity-aligned space for identity affine", {
  out <- output_aligned_space(list(c(5L, 6L, 7L), diag(4)))

  expect_equal(out$shape, c(5L, 6L, 7L))
  expect_equal(out$affine, diag(4))
  expect_equal(out$bounds$min, c(0, 0, 0))
  expect_equal(out$bounds$max, c(4, 5, 6))
})

test_that("output_aligned_space supports NeuroSpace and custom voxel sizes", {
  sp <- NeuroSpace(c(5L, 6L, 7L), spacing = c(2, 3, 4), origin = c(10, 20, 30))
  out <- output_aligned_space(sp, voxel_sizes = c(2, 3, 4))

  expect_equal(out$shape, c(5L, 6L, 7L))
  expect_equal(diag(out$affine)[1:3], c(2, 3, 4))
  expect_equal(out$affine[1:3, 4], c(10, 20, 30))
})

test_that("output_aligned_space handles >3D input by using first 3 dims", {
  expect_warning(
    out <- output_aligned_space(list(c(5L, 6L, 7L, 8L), diag(4))),
    "first three dimensions"
  )
  expect_equal(out$shape, c(5L, 6L, 7L))
})

test_that("vox2out_vox is an alias for output_aligned_space", {
  a <- output_aligned_space(list(c(4L, 4L, 4L), diag(4)))
  b <- vox2out_vox(list(c(4L, 4L, 4L), diag(4)))

  expect_equal(a$shape, b$shape)
  expect_equal(a$affine, b$affine)
})

test_that("slice_to_volume_affine supports R indexing", {
  aff <- slice_to_volume_affine(index = 3, axis = 3, shape = c(10, 10, 10), index_base = "R")
  expected <- matrix(
    c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 2,
      0, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  expect_equal(aff, expected)
})

test_that("slice_to_volume_affine supports zero-based axis/index compatibility", {
  aff <- slice_to_volume_affine(index = 2, axis = 2, shape = c(10, 10, 10), index_base = "zero")
  # axis=2 (zero-based z) => axis=3 in R, index0=2
  expect_equal(aff[3, 3], 2)
  expect_equal(dim(aff), c(4, 3))
})

test_that("slice2volume is an alias for slice_to_volume_affine", {
  a <- slice_to_volume_affine(index = 4, axis = 1, shape = c(10, 10, 10), index_base = "R")
  b <- slice2volume(index = 4, axis = 1, shape = c(10, 10, 10), index_base = "R")
  expect_equal(a, b)
})

test_that("slice_to_volume_affine validates bounds", {
  expect_error(slice_to_volume_affine(index = 0, axis = 1, shape = c(5, 5, 5), index_base = "R"))
  expect_error(slice_to_volume_affine(index = 5, axis = 0, shape = c(5, 5, 5), index_base = "zero"))
})
