test_that("SparseNeuroVec validity catches mismatched shapes", {
  sp <- NeuroSpace(c(8,8,8,5), c(2,2,2))
  bad_mask <- array(TRUE, c(4,4,4)) # too small relative to space
  dat <- matrix(0, nrow = 5, ncol = 10)
  expect_error(SparseNeuroVec(dat, sp, bad_mask), "Invalid mask dimensions|data argument must have three dimensions")
})

test_that("SparseNeuroVec validity enforces nrow(data)=time, ncol(data)=sum(mask)", {
  sp <- NeuroSpace(c(6,6,6,4), c(2,2,2))
  mask <- array(runif(6*6*6) > 0.8, c(6,6,6))
  dat_wrong_time <- matrix(0, nrow = 3, ncol = sum(mask))
  expect_error(SparseNeuroVec(dat_wrong_time, sp, mask), "Data/time mismatch")
  dat_wrong_cols <- matrix(0, nrow = 4, ncol = sum(mask) + 1)
  expect_error(SparseNeuroVec(dat_wrong_cols, sp, mask), "Data/mask mismatch|Matrix dimensions .* do not match mask cardinality")
})
