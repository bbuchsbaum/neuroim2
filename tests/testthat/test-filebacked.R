library(magrittr)
library(purrr)
library(testthat)
library(assertthat)

context("filebacked neurovec")

gmask <- system.file("extdata", "global_mask.nii", package="neuroim2")
gvec <- FileBackedNeuroVec("test_data/global_mask_v5.nii")

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

test_that("can concatenate two or more FileBackedNeuroVecs", {


})

test_that("can extract a single volume from a FileBackedNeuroVec", {
  vol1 <- drop(gvec[[1]])
  expect_equal(dim(gvec)[1:3], dim(vol1))
})

test_that("can extract a sub vector from a FileBackedNeuroVec", {
  vec1 <- sub_vector(gvec, 1:2)
  expect_true(inherits(vec1, "NeuroVec"))
  expect_equal(dim(vec1), c(dim(gvec)[1:3],2))
})

test_that("can map over each volume in a FileBackedNeuroVec", {
  mean.vol1 <- lapply(vols(gvec), mean)
  expect_equal(length(mean.vol1), dim(gvec)[4])
})

test_that("can map over first 50 vectors in a FileBackedNeuroVec", {
  mean.vec1 <- map_dbl(vectors(gvec, 1:50), mean)
  expect_equal(length(mean.vec1), 50)
})

test_that("can split a FileBackedNeuroVec into a set of clustered ROIs", {
  grid <- index_to_coord(space(gvec), as.numeric(1:prod(dim(gvec)[1:3])))
  kres <- kmeans(grid, centers=50)
  res <- gvec %>% split_clusters(kres$cluster) %>% map_dbl(~ mean(.))
  expect_equal(length(res), 50)

})

test_that("can map over a subset of vols in a FileBackedNeuroVec", {
  mean.vol1 <- lapply(vols(gvec, 2:3), mean)
  expect_equal(length(mean.vol1), 2)
})


test_that("can extract multiple series from a FileBackedNeuroVec", {
  expect_equal(as.vector(series(gvec, 1,1,1)), gvec[1,1,1,])
  mat <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))

  r1 <- apply(mat, 1, function(i) { series(gvec, i[1], i[2], i[3]) })
  r2 <- series(gvec, mat)
  expect_equal(r1, r2)
})

test_that("can extract an ROIVec from a FileBackedNeuroVec", {
  roi <- series_roi(gvec, 1:10)
  roi2 <- series_roi(gvec, LogicalNeuroVol(data=rep(1, 10), space=drop_dim(space(gvec)), indices=1:10))
  rvol <- ROIVol(space(drop_dim(space(gvec))), coords(roi2), data=rep(1, nrow(coords(roi2))))
  roi3 <- series_roi(gvec, rvol)
  roi4 <- series_roi(gvec, coords(rvol))
  expect_equal(coords(roi), coords(roi2))
  expect_equal(coords(roi2), coords(roi3))
  expect_equal(coords(roi3), coords(roi4))
})


test_that("can convert FileBackedNeuroVec to matrix", {

  mat <- as.matrix(gvec)

  ind <- sample(seq(1, prod(dim(gvec)[1:3])),100)
  mat2 <- t(do.call(cbind, lapply(ind, function(i) series(gvec, i))))

  expect_equal(mat[ind,], mat2)
})






