test_that("as.matrix(SparseNeuroVec) matches dense roundtrip", {
  set.seed(1)
  sp <- NeuroSpace(c(10,10,10,3), c(2,2,2))
  arr <- array(rnorm(prod(c(10,10,10,3))), dim = c(10,10,10,3))
  dense <- DenseNeuroVec(arr, sp)
  # Make a mask with ~30% voxels
  m <- array(runif(10*10*10) > 0.7, c(10,10,10))
  sparse <- as.sparse(dense, LogicalNeuroVol(m, drop_dim(space(dense))))
  md <- as.matrix(dense)
  ms <- as.matrix(sparse)
  # zero-out voxels outside mask in the dense matrix for fair comparison
  mask_vec <- as.vector(m)
  md[!mask_vec, ] <- 0
  expect_equal(dim(md), dim(ms))
  expect_equal(md, ms, tolerance = 1e-7)
})
