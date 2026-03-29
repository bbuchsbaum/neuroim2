test_that("bilateral_filter_4d_cpp_par is identity for zero window", {
  set.seed(1)
  arr <- array(rnorm(3 * 4 * 2 * 5), dim = c(3L, 4L, 2L, 5L))
  mask_idx <- which(array(TRUE, dim = dim(arr)[1:3]))

  out <- neuroim2:::bilateral_filter_4d_cpp_par(
    arr = arr,
    mask_idx = as.integer(mask_idx),
    spatial_window = 0L,
    temporal_window = 0L,
    spatial_sigma = 1,
    intensity_sigma = 1,
    temporal_sigma = 1,
    spacing = c(1, 1, 1, 1)
  )

  expect_equal(as.numeric(out), as.numeric(arr))
})

test_that("bilateral_filter_4d_cpp_par preserves constant arrays without NaNs", {
  arr <- array(5, dim = c(3L, 3L, 3L, 4L))
  mask_idx <- which(array(TRUE, dim = dim(arr)[1:3]))

  out <- neuroim2:::bilateral_filter_4d_cpp_par(
    arr = arr,
    mask_idx = as.integer(mask_idx),
    spatial_window = 1L,
    temporal_window = 1L,
    spatial_sigma = 1,
    intensity_sigma = 1,
    temporal_sigma = 1,
    spacing = c(1, 1, 1, 1)
  )

  expect_true(all(is.finite(as.numeric(out))))
  expect_equal(as.numeric(out), as.numeric(arr))
})

test_that("as.array works for SparseNeuroVec and preserves masked layout", {
  set.seed(2)
  d <- c(4L, 3L, 2L, 3L)
  sp <- NeuroSpace(d, spacing = c(1, 1, 1))
  arr <- array(rnorm(prod(d)), dim = d)
  dense <- DenseNeuroVec(arr, sp)

  mask_arr <- array(runif(prod(d[1:3])) > 0.4, dim = d[1:3])
  mask_arr[1] <- TRUE
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp))
  sparse <- as.sparse(dense, mask)

  sparse_arr <- as.array(sparse)
  expect_equal(dim(sparse_arr), d)

  expected_mat <- matrix(arr, nrow = prod(d[1:3]), ncol = d[4])
  expected_mat[!as.vector(mask_arr), ] <- 0
  expected_arr <- array(expected_mat, dim = d)

  expect_equal(sparse_arr, expected_arr, tolerance = 1e-8)
})

test_that("bilateral_filter_4d accepts SparseNeuroVec inputs", {
  set.seed(3)
  d <- c(4L, 4L, 3L, 4L)
  sp <- NeuroSpace(d, spacing = c(1, 1, 1))
  arr <- array(rnorm(prod(d)), dim = d)
  dense <- DenseNeuroVec(arr, sp)

  mask_arr <- array(runif(prod(d[1:3])) > 0.35, dim = d[1:3])
  mask_arr[1] <- TRUE
  mask <- LogicalNeuroVol(mask_arr, drop_dim(sp))
  sparse <- as.sparse(dense, mask)
  dense_from_sparse <- DenseNeuroVec(as.array(sparse), sp)

  args <- list(
    mask = mask,
    spatial_sigma = 1.2,
    intensity_sigma = 0.9,
    temporal_sigma = 1.0,
    spatial_window = 1L,
    temporal_window = 1L,
    temporal_spacing = 1
  )

  out_sparse <- do.call(bilateral_filter_4d, c(list(vec = sparse), args))
  out_dense <- do.call(bilateral_filter_4d, c(list(vec = dense_from_sparse), args))

  expect_s4_class(out_sparse, "DenseNeuroVec")
  expect_true(all(is.finite(as.numeric(as.array(out_sparse)))))
  expect_equal(as.array(out_sparse), as.array(out_dense), tolerance = 1e-8)
})
