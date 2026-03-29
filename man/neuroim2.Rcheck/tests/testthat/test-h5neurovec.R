library(purrr)

context("h5neurovec")

gmask <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")


gen_dat <- function(d1 = 12,
                    d2 = 12,
                    d3 = 12,
                    d4 = 4,
                    rand = FALSE) {
  if (rand) {
    dat <- array(rnorm(d1 * d2 * d3 * d4), c(d1, d2, d3, d4))
  } else {
    dat <- array(0, c(d1, d2, d3, d4))
  }
  spc <- NeuroSpace(c(d1, d2, d3, d4))
  DenseNeuroVec(dat, spc)
}


# test_that("can construct an H5NeuroVec", {
# 	bv <- array(rnorm(12*12*12*4), c(12,12,12,4))
# 	bv <- NeuroVec(bv, NeuroSpace(dim=c(12,12,12,4)))
# 	fname <- paste0(tempfile(), ".h5")
# 	hfile <- to_nih5_vec(bv,fname)
# 	hfile$close_all()
# 	hvec <- H5NeuroVec(fname)
#   expect_equal(hvec[1,1,1,1], bv[1,1,1,1],tolerance=1e-6)
#   expect_equal(hvec[1,1,1,], bv[1,1,1,],tolerance=1e-3)
#   expect_equal(hvec[1:5,1,1,1], bv[1:5,1,1,1],tolerance=1e-3)
#
# })


