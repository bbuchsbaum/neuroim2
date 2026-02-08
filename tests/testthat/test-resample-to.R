test_that("resample_to maps methods to interpolation integers", {
  skip_on_cran()
  img <- NeuroVol(array(rnorm(32*32*32), c(32,32,32)), NeuroSpace(c(32,32,32)))
  sp  <- NeuroSpace(c(24,24,24))
  expect_s4_class(resample_to(img, sp, method="nearest"), "NeuroVol")
  expect_s4_class(resample_to(img, sp, method="linear"),  "NeuroVol")
  expect_s4_class(resample_to(img, sp, method="cubic"),   "NeuroVol")
})

test_that("resample_to refuses unknown engine", {
  img <- NeuroVol(array(rnorm(16*16*16), c(16,16,16)), NeuroSpace(c(16,16,16)))
  sp  <- NeuroSpace(c(8,8,8))
  expect_error(resample_to(img, sp, engine="RNiftyReg"), "Only engine = 'internal'")
})

