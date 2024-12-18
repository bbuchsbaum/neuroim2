context("image processing")


test_that("can gaussian blur an image volume", {
  gmask <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
  gmask <- read_vol(gmask)
  g1 <- neuroim2::gaussian_blur(gmask, gmask, sigma=2)
  g2 <- neuroim2::gaussian_blur(gmask, gmask, sigma=8)
  g3 <- neuroim2::gaussian_blur(gmask, gmask, sigma=8, window=3)

  expect_true(!is.null(g1))
  expect_true(!is.null(g2))
  expect_true(!is.null(g3))

})

test_that("can apply guided filter image volume", {
  gmask <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")
  gmask <- read_vol(gmask)
  g1 <- neuroim2::guided_filter(gmask, 4)

  expect_true(!is.null(g1))

})
