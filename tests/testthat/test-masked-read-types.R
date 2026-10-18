library(testthat)
library(neuroim2)

# Masked reads of uncompressed and gzipped UINT/LONG/ULONG images must not
# wrap values above 2^31 (readBin(size = 4/8) used to sign-wrap them), whether
# or not the optional 'mmap' package is available.

write_typed <- function(x, datatype, ext) {
  f <- tempfile(fileext = ext)
  RNifti::writeNifti(x, f, datatype = datatype)
  f
}

check_masked <- function(x, datatype, ext) {
  f <- write_typed(x, datatype, ext)
  on.exit(unlink(f), add = TRUE)
  dense <- read_vec(f)
  expect_equal(as.vector(as.array(dense)), as.vector(x), info = paste(datatype, ext))
  m <- array(TRUE, dim(x)[1:3])
  m[1] <- FALSE
  mask <- LogicalNeuroVol(m, drop_dim(space(dense)))
  want <- matrix(x, ncol = dim(x)[4])[as.vector(m), , drop = FALSE]
  expect_warning(s <- read_vec(f, mask = mask), NA)
  expect_s4_class(s, "SparseNeuroVec")
  expect_equal(unname(series(s, which(m))), unname(t(want)), info = paste(datatype, ext))
  # A subset of volumes, out of order and repeated.
  s2 <- read_vec(f, mask = mask, indices = c(2L, 1L, 2L))
  expect_equal(unname(series(s2, which(m))), unname(t(want[, c(2, 1, 2)])),
               info = paste(datatype, ext, "indices"))
}

d <- c(3, 3, 2, 2)

test_that("masked uint32 reads keep values above 2^31", {
  set.seed(2)
  x <- array(round(runif(prod(d), 0, 4e9)), d)
  x[2] <- 2809496144
  for (ext in c(".nii", ".nii.gz")) check_masked(x, "uint32", ext)
})

test_that("masked int64 reads keep values beyond the 32-bit range", {
  x <- array(seq(-5, 30, length.out = prod(d)), d)
  x[2] <- 3e9
  x[3] <- -3e9
  for (ext in c(".nii", ".nii.gz")) check_masked(x, "int64", ext)
})

test_that("masked wide-integer reads agree when mmap is unavailable", {
  x <- array(round(seq(1, 4e9, length.out = prod(d))), d)
  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  check_masked(x, "uint32", ".nii")
})
