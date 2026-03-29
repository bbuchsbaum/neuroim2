context("io-scaling-and-truncation")

test_that("read_elements() errors on unexpected EOF", {
  tmp <- tempfile(fileext = ".bin")
  on.exit(unlink(tmp))

  con <- file(tmp, open = "wb")
  writeBin(as.double(1:5), con, size = 8)
  close(con)

  reader <- BinaryReader(tmp, byte_offset = 0L, data_type = "double", bytes_per_element = 8L)
  on.exit(close(reader), add = TRUE)

  expect_error(read_elements(reader, 10), "Unexpected EOF")
})

test_that("NeuroVolSource applies slope/intercept scaling", {
  f3d <- tempfile(fileext = ".nii")
  on.exit(unlink(f3d), add = TRUE)

  arr3d <- array(seq_len(8), dim = c(2, 2, 2))
  RNifti::writeNifti(arr3d, f3d)

  meta <- read_header(f3d)
  meta@slope <- 2
  meta@intercept <- 10

  src <- new("NeuroVolSource", meta_info = meta, index = 1L)
  vol <- load_data(src)

  expect_equal(as(vol, "array"), array(as.numeric(arr3d) * 2 + 10, dim = dim(arr3d)))
})

test_that("NeuroVecSource applies per-volume slope/intercept scaling", {
  f4d <- tempfile(fileext = ".nii")
  on.exit(unlink(f4d), add = TRUE)

  arr4d <- array(seq_len(24), dim = c(2, 2, 2, 3))
  RNifti::writeNifti(arr4d, f4d)

  meta <- read_header(f4d)
  meta@slope <- c(1, 2, 3)
  meta@intercept <- c(0, 10, 100)

  src <- new("NeuroVecSource", meta_info = meta, indices = as.integer(c(1, 3)))
  vec <- load_data(src)

  out <- as(vec, "array")
  expect_equal(dim(out), c(2, 2, 2, 2))

  expect_equal(out[, , , 1], array(as.numeric(arr4d[, , , 1]) * 1 + 0, dim = c(2, 2, 2)))
  expect_equal(out[, , , 2], array(as.numeric(arr4d[, , , 3]) * 3 + 100, dim = c(2, 2, 2)))
})

test_that("SparseNeuroVecSource reads only masked voxels and applies scaling", {
  f4d <- tempfile(fileext = ".nii")
  on.exit(unlink(f4d), add = TRUE)

  arr4d <- array(seq_len(24), dim = c(2, 2, 2, 3))
  RNifti::writeNifti(arr4d, f4d)

  meta <- read_header(f4d)
  meta@slope <- c(1, 2, 3)
  meta@intercept <- c(0, 10, 100)

  mask <- array(FALSE, dim = c(2, 2, 2))
  mask[1, 1, 1] <- TRUE
  mask[2, 2, 2] <- TRUE
  src <- neuroim2:::SparseNeuroVecSource(meta, indices = as.integer(c(2, 3)), mask = mask)
  svec <- load_data(src)

  idx <- indices(svec)

  expect_equal(nrow(svec@data), 2)
  expect_equal(ncol(svec@data), length(idx))

  exp2 <- as.numeric(as.vector(arr4d[, , , 2])[idx]) * 2 + 10
  exp3 <- as.numeric(as.vector(arr4d[, , , 3])[idx]) * 3 + 100
  expect_equal(svec@data[1, ], exp2)
  expect_equal(svec@data[2, ], exp3)
})
