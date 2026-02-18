context("hypervec io")

test_that("read_hyper_vec reads 5D nifti into NeuroHyperVec", {
  arr5 <- array(seq_len(2 * 2 * 2 * 3 * 4), dim = c(2, 2, 2, 3, 4))
  f5d <- tempfile(fileext = ".nii")
  on.exit(unlink(f5d), add = TRUE)

  RNifti::writeNifti(arr5, f5d)

  hvec <- read_hyper_vec(f5d)
  expect_s4_class(hvec, "NeuroHyperVec")
  expect_equal(dim(hvec), dim(arr5))
  expect_equal(neuroim2:::dense_array_5d(hvec), arr5)
})

test_that("read_hyper_vec applies optional spatial mask", {
  arr5 <- array(seq_len(2 * 2 * 2 * 2 * 2), dim = c(2, 2, 2, 2, 2))
  f5d <- tempfile(fileext = ".nii")
  on.exit(unlink(f5d), add = TRUE)

  RNifti::writeNifti(arr5, f5d)

  mask <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 2))
  hvec <- read_hyper_vec(f5d, mask = mask)
  dense <- neuroim2:::dense_array_5d(hvec)

  expected <- arr5
  for (m in seq_len(dim(arr5)[5])) {
    for (l in seq_len(dim(arr5)[4])) {
      vol <- expected[, , , l, m]
      vol[!mask] <- 0
      expected[, , , l, m] <- vol
    }
  }

  expect_equal(dense, expected)
  expect_equal(sum(hvec@mask@.Data), sum(mask))
})

test_that("read_image auto-dispatches 5D inputs to NeuroHyperVec", {
  arr5 <- array(seq_len(2 * 2 * 2 * 2 * 3), dim = c(2, 2, 2, 2, 3))
  f5d <- tempfile(fileext = ".nii")
  on.exit(unlink(f5d), add = TRUE)

  RNifti::writeNifti(arr5, f5d)

  img <- read_image(f5d)
  expect_s4_class(img, "NeuroHyperVec")
  expect_equal(dim(img), dim(arr5))
})

test_that("write_vec writes NeuroHyperVec to 5D nifti", {
  arr5 <- array(seq_len(2 * 2 * 2 * 2 * 3), dim = c(2, 2, 2, 2, 3))
  fin <- tempfile(fileext = ".nii")
  fout <- tempfile(fileext = ".nii")
  on.exit(unlink(c(fin, fout)), add = TRUE)

  RNifti::writeNifti(arr5, fin)
  hvec <- read_hyper_vec(fin)

  write_vec(hvec, fout)

  hdr <- read_header(fout)
  expect_equal(dim(hdr), dim(arr5))

  arr_out <- as.array(RNifti::readNifti(fout))
  if (!identical(dim(arr_out), dim(arr5))) {
    arr_out <- array(as.numeric(arr_out), dim = dim(arr5))
  }
  expect_equal(dim(arr_out), dim(arr5))
  expect_equal(as.numeric(arr_out), as.numeric(arr5))
})
