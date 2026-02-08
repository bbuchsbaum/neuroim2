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
