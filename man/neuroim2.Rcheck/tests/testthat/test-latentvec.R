# library(purrr)
# gmask <- system.file("extdata", "global_mask.nii", package="neuroim2")
#
#
# gen_dat <- function(d1 = 12,
#                     d2 = 12,
#                     d3 = 12,
#                     d4 = 4,
#                     rand = FALSE) {
#   if (rand) {
#     dat <- array(rnorm(d1 * d2 * d3 * d4), c(d1, d2, d3, d4))
#   } else {
#     dat <- array(0, c(d1, d2, d3, d4))
#   }
#   spc <- NeuroSpace(c(d1, d2, d3, d4))
#   DenseNeuroVec(dat, spc)
# }
#
# gen_latent_vec <- function(d1,d2,d3,d4) {
#   bv <- gen_dat(12,12,12,4, rand=TRUE)
#   mat <- bv@.Data
#   dim(mat) <- c(12*12*12,4)
#   mat <- t(mat)
#   mask <- bv[[1]]
#   mask[] <- 1
#   mask <- as.logical(mask)
#   pres <- prcomp(mat)
#   svec <- LatentNeuroVec(pres$x, pres$rotation, space=space(bv), mask=mask, offset=colMeans(mat))
#
# }
#
# context("latentneurovec")
#
# test_that("can construct a LatentNeuroVec", {
#   bv <- gen_dat(12,12,12,4, rand=TRUE)
#   mat <- bv@.Data
#   dim(mat) <- c(12*12*12,4)
#   mat <- t(mat)
#   mask <- bv[[1]]
#   mask[] <- 1
#   mask <- as.logical(mask)
#   pres <- prcomp(mat)
#   svec <- LatentNeuroVec(pres$x, pres$rotation, space=space(bv), mask=mask, offset=colMeans(mat))
#
#   expect_true(!is.null(svec))
#   expect_equal(svec[1,1,1,1], bv[1,1,1,1])
#   expect_equal(bv[12,12,12,4], svec[12,12,12,4])
#   expect_equal(dim(svec), c(12,12,12,4))
#   expect_equal(series(svec, 1:2), series(bv,1:2))
#
#   ## TODO should series drop dimensions?
#   expect_equal(as.vector(series(svec, 3)), series(bv,3))
#
# })
#
# # test_that("can write a LatentNeuroVec to h5", {
# #   bv <- gen_dat(12,12,12,4, rand=TRUE)
# #   mat <- bv@.Data
# #   dim(mat) <- c(12*12*12,4)
# #   mat <- t(mat)
# #   mask <- bv[[1]]
# #   mask[] <- 1
# #   mask <- as.logical(mask)
# #   pres <- prcomp(mat)
# #   svec <- LatentNeuroVec(pres$x, pres$rotation, space=space(bv), mask=mask, offset=colMeans(mat))
# #   tmp <- paste0(tempfile())
# #   write_vec(svec, tmp)
# # })
#
# test_that("can extract a single volume from a LatentNeuroVec", {
#   bv1 <- gen_latent_vec()
#   bv2 <- gen_latent_vec()
#
#   bv3 <- concat(bv1, bv2)
#
#   expect_equal(as.vector(series(bv1,1)), series(bv3,1)[1:4])
#   expect_equal(as.vector(series(bv2,1)), series(bv3,1)[5:8])
# })
