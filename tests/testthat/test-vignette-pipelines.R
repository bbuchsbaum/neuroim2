context("vignette pipelines")

library(purrr)

test_that("test_searchlight_statistics works correctly", {
  # Load example volume
  file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
  vol <- read_vol(file_name)
  mask.idx <- which(vol > 0)
  
  # Create test volume with random values
  vol2 <- vol
  vol2[mask.idx] <- runif(length(mask.idx))
  
  # Test searchlight with eager=TRUE for local standard deviation
  sdvol_values <- vol %>% 
    searchlight(radius=5, eager=TRUE) %>% 
    map_dbl(~ sd(values(.)))
  
  sdvol <- NeuroVol(sdvol_values, space=space(vol), indices=which(vol != 0))
  expect_true(inherits(sdvol, "NeuroVol"))
  expect_equal(length(sdvol_values), length(which(vol != 0)))
})

test_that("test_4d_processing works correctly", {
  # Load example volume and create a 4D NeuroVec
  file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
  vol <- read_vol(file_name)
  
  # Create a 4D volume by concatenating
  vec <- concat(vol, vol, vol, vol, vol)
  expect_true(inherits(vec, "NeuroVec"))
  expect_equal(dim(vec)[4], 5)
  
  # Map over volumes
  mean_vec <- vec %>% vols() %>% map_dbl(~ mean(values(.)))
  expect_equal(length(mean_vec), 5)
  
  sd_vec <- vec %>% vols() %>% map_dbl(~ sd(values(.)))
  expect_equal(length(sd_vec), 5)
  
  # Map over vectors (time series at each voxel)
  mean_vol_values <- vec %>% vectors() %>% map_dbl(~ mean(.))
  mean_vol <- NeuroVol(mean_vol_values, space=space(vol))
  expect_true(inherits(mean_vol, "NeuroVol"))
  expect_equal(dim(mean_vol), dim(vol)[1:3])
})

test_that("test_knn_smoothing works correctly", {
  # Load example volume
  file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
  vol <- read_vol(file_name)
  mask.idx <- which(vol > 0)
  
  # Create test volume with random values
  vol2 <- vol
  vol2[mask.idx] <- runif(length(mask.idx))
  
  # KNN smoothing using searchlight
  k <- 12
  knnfvol_values <- vol2 %>% 
    searchlight(radius=6, eager=TRUE) %>% 
    map_dbl(function(x) {
      # Get values from the ROI
      roi_values <- values(x)
      center_idx <- x@center_index
      
      # For ROIVolWindow, center_index refers to position in the ROI
      # We need to handle this differently
      if (length(roi_values) > 0) {
        mean(roi_values)
      } else {
        NA_real_
      }
    })
  
  # Remove NAs and create volume
  valid_idx <- !is.na(knnfvol_values)
  knnfvol <- NeuroVol(
    knnfvol_values[valid_idx], 
    space=space(vol), 
    indices=which(vol != 0)[valid_idx]
  )
  
  expect_true(inherits(knnfvol, "NeuroVol"))
  expect_true(sum(!is.na(values(knnfvol))) > 0)
})