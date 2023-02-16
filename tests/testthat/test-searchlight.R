context("searchlight")

gmask <- system.file("extdata", "global_mask.nii", package="neuroim2")
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
  expect_equal(length(s1), 200)
})


