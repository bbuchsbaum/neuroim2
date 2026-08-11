test_that("as.raster works for NeuroSlice", {
  space <- NeuroSpace(c(10, 10))
  dat <- matrix(rnorm(100), 10, 10)
  sl <- NeuroSlice(dat, space)
  r <- as.raster(sl)
  expect_s3_class(r, "raster")
  expect_equal(dim(r), dim(dat))
})

test_that("as.raster works for NeuroVol", {
  space <- NeuroSpace(c(10,10,5))
  dat <- array(rnorm(10*10*5), c(10,10,5))
  vol <- NeuroVol(dat, space)
  r <- as.raster(vol, zlevel=1)
  expect_s3_class(r, "raster")
  expect_equal(dim(r), dim(dat)[1:2])
})
