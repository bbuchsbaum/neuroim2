context("searchlight")

gmask <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
gmask_gz <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")


test_that("can extract searchlight object from a 3d volume", {
  vol <- read_vol(gmask)
  s1 <- searchlight(vol, 8, eager=FALSE, nonzero=FALSE)
  s2 <- searchlight(vol, 8, eager=TRUE,  nonzero=FALSE)

  s3 <- searchlight(vol, 8, eager=FALSE,  nonzero=TRUE)
  s4 <- searchlight(vol, 8, eager=TRUE,  nonzero=TRUE)

  expect_equal(nrow(coords(s1[[1]])), nrow(coords(s2[[1]])))
  expect_equal(nrow(coords(s3[[1]])), nrow(coords(s4[[1]])))
})

test_that("can extract clustered searchlight object from a 3d volume", {
  vol <- read_vol(gmask)
  s1 <- clustered_searchlight(vol, csize=100)
  expect_equal(length(s1), 100)

  grid <- index_to_grid(vol, which(vol>0))
  kres <- kmeans(grid, centers=200, iter.max=500)
  cvol <- ClusteredNeuroVol(vol, clusters=kres$cluster)
  s2 <- clustered_searchlight(vol, cvol=cvol)
  expect_equal(length(s2), 200)
})


test_that("random_searchlight works as expected", {
  mask <- read_vol(gmask)
  random_sl <- random_searchlight(mask, radius = 5)
  expect_true(!is.null(random_sl), info = "random_searchlight failed to create an iterator")
})

test_that("bootstrap_searchlight works as expected", {
  mask <- read_vol(gmask)
  bootstrap_sl <- bootstrap_searchlight(mask, radius = 8, iter = 100)
  expect_true(!is.null(bootstrap_sl), info = "bootstrap_searchlight failed to create an iterator")
})

test_that("searchlight_coords works as expected", {
  mask <- read_vol(gmask)
  coords_list <- searchlight_coords(mask, radius = 4, nonzero = FALSE, cores = 0)
  expect_true(!is.null(coords_list), info = "searchlight_coords failed to create an iterator")
})

test_that("searchlight works as expected", {
  mask <- read_vol(gmask)
  searchlight_list <- searchlight(mask, radius = 4, eager = FALSE, nonzero = FALSE, cores = 0)
  expect_true(!is.null(searchlight_list), info = "searchlight failed to create an iterator")
})

test_that("clustered_searchlight works as expected", {
  mask <- read_vol(gmask)
  clustered_sl <- clustered_searchlight(mask, csize = 5)
  expect_true(!is.null(clustered_sl), info = "clustered_searchlight failed to create an iterator")
})



