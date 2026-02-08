context("vignette neuro vectors")

test_that("test_create_neurovec works correctly", {
  skip_on_cran()
  # Load example data
  file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
  
  # Read single file as NeuroVec
  vec <- read_vec(file_name)
  expect_true(inherits(vec, "NeuroVec"))
  expect_equal(length(dim(vec)), 4)
  
  # Read multiple files
  vec_multi <- read_vec(c(file_name, file_name, file_name))
  expect_equal(dim(vec_multi)[4], 12)  # 3 files * 4 volumes each
  
  # Test sub_vector
  vec_subset <- sub_vector(vec_multi, 1:6)
  expect_equal(dim(vec_subset)[4], 6)
  
  # Test series extraction
  ts <- series(vec, 1, 1, 1)
  expect_true(is.numeric(ts))
  expect_equal(length(ts), dim(vec)[4])
  
  # Test series_roi with spherical ROI
  vol <- read_vol(file_name)
  roi <- spherical_roi(vol, c(12, 12, 12), radius=8)
  rvec <- series_roi(vec, roi)
  expect_true(inherits(rvec, "ROIVec"))
  
  # Test series_roi with indices
  rvec_idx <- series_roi(vec, 1:100)
  expect_true(inherits(rvec_idx, "ROIVec"))
  expect_equal(nrow(coords(rvec_idx)), 100)
  
  # Test series with indices
  mat <- series(vec, 1:100)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), dim(vec)[4])  # series returns time x voxels
  expect_equal(ncol(mat), 100)
  
  # Test series with voxel coordinates
  mask <- read_vol(file_name)
  vox <- index_to_grid(mask, 1:100)
  mat_vox <- series(vec, vox)
  expect_true(is.matrix(mat_vox))
  expect_equal(nrow(mat_vox), dim(vec)[4])  # series returns time x voxels
  expect_equal(ncol(mat_vox), 100)
})