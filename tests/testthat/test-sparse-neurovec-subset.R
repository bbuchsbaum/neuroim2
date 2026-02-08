context("sparse neurovec subset")

library(neuroim2)

make_sparse_vec2 <- function() {
  sp <- NeuroSpace(c(2,2,2,2))
  mask_arr <- array(c(1,0,1,0, 0,1,0,1), dim=c(2,2,2))
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,2)))
  nvox <- sum(mask_arr)
  # data: time x voxels
  data <- matrix(seq_len(2 * nvox), nrow = 2, ncol = nvox)
  SparseNeuroVec(data, sp, mask)
}

test_that("[ on AbstractSparseNeuroVec matches dense values and drop behavior", {
  svec <- make_sparse_vec2()
  dvec <- as(svec, "DenseNeuroVec")

  # multi-voxel subset retains shape
  sub_s <- svec[1:2, 1:2, 1:2, 1:2, drop = FALSE]
  sub_d <- dvec[1:2, 1:2, 1:2, 1:2, drop = FALSE]
  expect_equal(sub_s, sub_d)

  # single voxel with drop
  s_drop <- svec[1,1,1,1, drop = TRUE]
  d_drop <- dvec[1,1,1,1, drop = TRUE]
  expect_equal(s_drop, d_drop)

  # single voxel without drop keeps length-one time dim
  s_nodrop <- svec[1,1,1,1, drop = FALSE]
  d_nodrop <- dvec[1,1,1,1, drop = FALSE]
  expect_equal(dim(s_nodrop), c(1,1,1,1))
  expect_equal(s_nodrop, d_nodrop)
})
