context("sparse neurovec regressions")

library(neuroim2)

make_sparse_vec <- function() {
  sp <- NeuroSpace(c(2,2,2,2))
  mask_arr <- array(c(1,0,1,0, 0,1,0,1), dim=c(2,2,2))
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,2)))
  nvox <- sum(mask_arr)
  data <- matrix(seq_len(2 * nvox), nrow = 2, ncol = nvox)  # time x voxels
  SparseNeuroVec(data, sp, mask)
}

test_that("[ subsetting matches dense conversion", {
  svec <- make_sparse_vec()
  dvec <- as(svec, "DenseNeuroVec")
  expect_equal(svec[1:2,1:2,1:2,1:2], dvec[1:2,1:2,1:2,1:2])
  expect_equal(svec[1,1,1,1, drop=TRUE], dvec[1,1,1,1, drop=TRUE])
})

test_that("linear_access returns zeros for masked-out elements", {
  svec <- make_sparse_vec()
  # choose a masked-out spatial index: (2,2,1) -> linear 4
  idx <- 4  # in 3D order
  # timepoint 1 => overall linear index: spatial 4 at time1
  lin <- idx + (1-1)*prod(dim(svec)[1:3])
  expect_equal(linear_access(svec, lin), 0)
})

test_that("as.dense populates only masked voxels", {
  svec <- make_sparse_vec()
  dvec <- as(svec, "DenseNeuroVec")
  mask_idx <- which(svec@mask@.Data)
  m <- matrix(dvec@.Data, nrow = prod(dim(dvec)[1:3]), ncol = dim(dvec)[4])
  expect_true(all(m[-mask_idx, ] == 0))
  expect_true(any(m[mask_idx, ] != 0))
})
