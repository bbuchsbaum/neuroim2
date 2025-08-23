context("vignette roi")

library(purrr)

test_that("test_searchlight works correctly", {
  # Load example volume
  file_name <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
  vol <- read_vol(file_name)
  
  # Test exhaustive searchlight with eager=TRUE
  slist <- searchlight(vol, radius=8, eager=TRUE)
  expect_true(is.list(slist))
  expect_true(length(slist) > 0)
  
  # Compute mean values in each searchlight
  ret <- slist %>% purrr::map_dbl(~ mean(vol[coords(.)]))
  expect_equal(length(ret), length(slist))
  expect_true(all(!is.na(ret)))
  
  # Test random searchlight
  ret_random <- vol %>% 
    random_searchlight(radius=8) %>% 
    purrr::map_dbl(~ mean(vol[coords(.)]))
  expect_true(length(ret_random) > 0)
  expect_true(all(!is.na(ret_random)))
  
  # Test clustered searchlight
  grid <- index_to_coord(vol, which(vol > 0))
  kres <- kmeans(grid, centers=50, iter.max=500)
  kvol <- ClusteredNeuroVol(vol, kres$cluster)
  
  ret_clustered <- vol %>% 
    clustered_searchlight(cvol=kvol) %>% 
    purrr::map_dbl(~ mean(vol[coords(.)]))
  expect_equal(length(ret_clustered), 50)
  expect_true(all(!is.na(ret_clustered)))
  
  # Test searchlight_coords
  coords_list <- searchlight_coords(vol, radius=4, nonzero=FALSE)
  expect_true(inherits(coords_list, "deflist"))
  
  # Force evaluation of a few elements to test
  first_coords <- coords_list[[1]]
  expect_true(is.matrix(first_coords))
  expect_equal(ncol(first_coords), 3)
})