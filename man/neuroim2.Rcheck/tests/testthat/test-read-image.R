context("read_image")

test_that("read_image returns NeuroVol for 3D files", {
  arr3d <- array(seq_len(8), dim = c(2, 2, 2))
  f3d <- tempfile(fileext = ".nii")
  on.exit(unlink(f3d))

  RNifti::writeNifti(arr3d, f3d)

  img <- read_image(f3d)
  expect_s4_class(img, "NeuroVol")
  expect_equal(dim(img), c(2, 2, 2))
  expect_equal(as(img, "array"), arr3d)
})

test_that("read_image returns NeuroVec for 4D files and respects overrides", {
  arr4d <- array(seq_len(24), dim = c(2, 2, 2, 3))
  f4d <- tempfile(fileext = ".nii")
  on.exit(unlink(f4d), add = TRUE)

  RNifti::writeNifti(arr4d, f4d)

  vec <- read_image(f4d)
  expect_s4_class(vec, "NeuroVec")
  expect_equal(dim(vec), c(2, 2, 2, 3))

  vec_slice <- read_image(f4d, index = 2)
  expect_equal(dim(vec_slice), c(2, 2, 2, 1))
  expect_equal(as(vec_slice, "array")[, , , 1], arr4d[, , , 2])

  vol <- read_image(f4d, type = "vol", index = 3)
  expect_s4_class(vol, "NeuroVol")
  expect_equal(dim(vol), c(2, 2, 2))
  expect_equal(as(vol, "array"), arr4d[, , , 3])
})
