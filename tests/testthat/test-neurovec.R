
library(purrr)
gmask <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")


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

context("neurovec")

test_that("can construct a DenseNeuroVec", {
	bv <- gen_dat(12,12,12,4)
	expect_true(!is.null(bv))
	expect_equal(bv[1,1,1,1], 0)
	expect_equal(bv[12,12,12,4], 0)
	expect_equal(dim(bv), c(12,12,12,4))
})

test_that("can construct a DenseNeuroVec from a matrix", {
  spc <- NeuroSpace(c(10,10,10,4))
  mat <- matrix(rnorm(4*10^3), 4, 10^3)
  vec1 <- NeuroVec(t(mat),spc)
  vec2 <- NeuroVec(mat,spc)
  expect_equal(vec1,vec2)
  expect_equal(length(vec1),length(vec2))
})

test_that("can construct a SparseNeuroVec from NeuroVec constructor", {
  mask <- read_vol(gmask)
  idx <- which(mask>0)
  mat <- matrix(rnorm(length(idx))*5, length(idx), 5)
  vec <- NeuroVec(mat,add_dim(space(mask),5), mask=as.logical(mask))
  expect_equal(dim(vec), c(dim(mask),5))
})

test_that("can construct a DenseNeuroVec from a list of vols", {
  spc <- NeuroSpace(c(10,10,10))
  vol <- NeuroVol(rnorm(10*10*10), spc)
  vlist <- list(vol,vol,vol)
  vec <- NeuroVec(vlist)
  expect_equal(dim(vec), c(10,10,10,3))
})

test_that("can concatenate two or more NeuroVecs", {
	bv1 <- gen_dat(rand=TRUE)
	bv2 <- gen_dat(rand=TRUE)

	bv3 <- concat(bv1, bv2)
	expect_true(inherits(bv3, "NeuroVec"))
	expect_equal(dim(bv3), c(12,12,12,8))

	bv4 <- concat(bv1,bv2, bv1, bv3)
	expect_equal(dim(bv4), c(12,12,12,20))
	expect_true(inherits(bv4, "NeuroVec"))
	expect_equal(bv4[1,1,1,1],bv1[1,1,1,1])

})

test_that("can extract a single volume from a NeuroVec", {
	bv1 <- gen_dat()
	bv2 <- gen_dat()

	bv3 <- concat(bv1, bv2)

	vol1 <- bv3[[1]]
	expect_equal(dim(vol1), c(12,12,12))
})


test_that("can extract a sub vector from a NeuroVec", {
  bv1 <- gen_dat()
  vec1 <- sub_vector(bv1, 1:2)
  expect_true(inherits(vec1, "NeuroVec"))
  expect_equal(dim(vec1), c(12,12,12,2))
})

test_that("can map over each volume in a NeuroVec", {
	bv1 <- gen_dat(rand=TRUE)
	mean.vol1 <- lapply(vols(bv1), mean)
	expect_equal(length(mean.vol1), 4)
	expect_equal(unlist(mean.vol1), apply(bv1, 4, mean))
})

test_that("can map over each vector in a NeuroVec", {
  bv1 <- gen_dat(10,10,10, rand=TRUE)
  mean.vec1 <- map_dbl(vectors(bv1), mean)
  mean.vec2 <- map_dbl(vectors(bv1, 1:10), mean)
  expect_equal(length(mean.vec1), 10*10*10)
  expect_equal(length(mean.vec2), 10)
  expect_equal(unlist(mean.vec1), as.vector(apply(bv1, 1:3, mean)))
})

test_that("can split a NeuroVec into a set of subvectors", {
  bv1 <- gen_dat(10,10,10, d4=20, rand=TRUE)
  blstruc <- rep(1:4, length.out=20)

  blks <- split_blocks(bv1, blstruc)
  assert_that(length(blks) == length(unique(blstruc)))
  assert_that(dim(blks[[1]])[4] == sum(blstruc == 1))

  ## nested map
  ## first divide into blocks, then convert to vectors, then foreach block/vector compute mean
  res <- bv1 %>% split_blocks(blstruc) %>% map(vectors) %>% map(~ map_dbl(., mean))
  expect_equal(length(res), length(blks))
  expect_equal(length(res[[1]]), prod(dim(bv1)[1:3]))

})

test_that("can split a NeuroVec into a set of clustered ROIs", {
  bv1 <- gen_dat(10,10,10, d4=20, rand=TRUE)
  blstruc <- rep(1:4, length.out=20)

  grid <- index_to_coord(space(bv1), as.numeric(1:prod(dim(bv1)[1:3])))
  kres <- kmeans(grid, centers=50)
  res <- bv1 %>% split_clusters(kres$cluster) %>% map_dbl(~ mean(.))
  expect_true(TRUE)
})


test_that("can extract a list of vectors from a NeuroVec", {
  bv1 <- gen_dat(10,10,10, rand=TRUE)

  vecs <- bv1 %>% vectors
  expect_equal(length(vecs), prod(dim(bv1)[1:3]))
  expect_equal(length(vecs[[1]]), dim(bv1)[4])

  vecs2 <- bv1 %>% vectors(subset=1:100)
  expect_equal(length(vecs2), 100)
  expect_equal(length(vecs2[[1]]), dim(bv1)[4])

  keep <- runif(prod(dim(bv1)[1:3])) > .5
  vecs3 <- bv1 %>% vectors(subset=keep)

  expect_equal(length(vecs3), sum(keep))
  expect_equal(length(vecs3[[1]]), dim(bv1)[4])

})


test_that("can map over a subset of vols in a NeuroVec", {
  bv1 <- gen_dat()
  mean.vol1 <- lapply(vols(bv1, 2:3), mean)
  expect_equal(length(mean.vol1), 2)
  expect_equal(unlist(mean.vol1), apply(bv1, 4, mean)[2:3])
})


test_that("can extract multiple series from a NeuroVec", {
	bv1 <- gen_dat()
	expect_equal(series(bv1, 1,1,1), bv1[1,1,1,])
  mat <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))

	r1 <- apply(mat, 1, function(i) { series(bv1, i[1], i[2], i[3]) })
	r2 <- series(bv1, mat)
	expect_equal(r1, r2)
})

test_that("can extract an ROIVec from a NeuroVec", {
  bv1 <- gen_dat()
  roi <- series_roi(bv1, 1:10)
  roi2 <- series_roi(bv1, LogicalNeuroVol(data=rep(1, 10), space=drop_dim(space(bv1)), indices=1:10))
  rvol <- ROIVol(space(drop_dim(space(bv1))), coords(roi2), data=rep(1, nrow(coords(roi2))))
  roi3 <- series_roi(bv1, rvol)
  roi4 <- series_roi(bv1, coords(rvol))
  expect_equal(coords(roi), coords(roi2))
  expect_equal(coords(roi2), coords(roi3))
  expect_equal(coords(roi3), coords(roi4))
})


test_that("can convert NeuroVec to matrix", {
	bv1 <- gen_dat(5,5,5,5, rand=TRUE)
	mat <- as(bv1, "matrix")

	ind <- 1:(5*5*5)
	mat2 <- t(do.call(cbind, lapply(ind, function(i) series(bv1, i))))

	expect_equal(mat, mat2)
})

test_that("can perform arithmetic on NeuroVec", {
  bv1 <- gen_dat(5,5,5,5)
  bv2 <- gen_dat(5,5,5,5)
  bv3 <- bv1+bv2
  bv4 <- bv1*bv2
  bv5 <- bv2-bv1
  expect_true(TRUE)
})



# test.SparseNeuroVec.series <- function() {
#
# 	spvec <- loadVector("data/qrscan01.nii.gz", indices=1:4, mask=rep(TRUE, 96*96*26))
# 	bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4)
#
# 	voxmat <- rbind(c(10,10,10), c(20,20,10), c(30,30, 10), c(40,40,10), c(50,50,10))
#
# 	expect_equal(series(bvec, voxmat), series(spvec, voxmat))
# }



# test.SparseNeuroVec.roundtrip.io <- function() {
# 	bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4, mask=rep(TRUE, 96*96*26))
# 	template <- takeVolume(bvec,1)
#
#
# 	mask.idx <- sort(sample(1:length(template), 1000))
# 	vals <- rnorm(length(mask.idx))
# 	bv <- NeuroVol(vals, space(template), indices=mask.idx)
# 	fname <- paste(tempfile(), ".nii", sep="")
# 	writeVolume(bv,fname)
# 	bv2 <- loadVolume(fname)
#
# 	expect_equal(dim(bv2), dim(bv))
# 	expect_equalNumeric(trans(bv2), trans(bv), tol=.0001)
#
#
# 	expect_equal(dim(bv2), dim(bv))
# 	expect_equalNumeric(trans(bv2), trans(bv), tol=.0001)
#
# }

test_that("can convert dense NeuroVec to sparse", {
	bv1 <- gen_dat()
	svec <- as.sparse(bv1, c(1,100,1000))
	expect_equal(dim(svec), c(12,12,12,4))

})



test_that("can construct a SparseNeuroVec", {
	dat <- array(rnorm(12*12*12*4), c(12,12,12,4))
	spc <- NeuroSpace(c(12,12,12,4))
	tmp <- rnorm(12*12*12)
	mask <- tmp > -1000000
	mask <- LogicalNeuroVol(mask, drop_dim(spc))
	bvec <- SparseNeuroVec(dat, spc, mask)

	expect_equal(dim(bvec), dim(dat))
	expect_equal(dat[1,1,1,1], bvec[1,1,1,1])
	expect_equal(dat[1,1,1,1], bvec[1,1,1,1])
	expect_equal(dat[1,1,1,], bvec[1,1,1,])
	expect_equal(dat[1,1,,], bvec[1,1,,])
	expect_equal(dat[1,,,], bvec[1,,,])
	expect_equal(dat[1,,,], bvec[1,,,])
	expect_equal(dat[1,2:3,,], bvec[1,2:3,,])
	expect_equal(dat[1,2:3,,2:3], bvec[1,2:3,,2:3])
	expect_equal(dat[1:3,2:3,,], bvec[1:3,2:3,,])
	expect_equal(dat[1,2:3,2:3,], bvec[1,2:3,2:3,])

})



test_that("can perform arithmetic on a SparseNeuroVec", {
  dat <- array(rnorm(12*12*12*4), c(12,12,12,4))
  spc <- NeuroSpace(c(12,12,12,4))
  tmp <- array(rnorm(12*12*12), c(12,12,12))
  mask <- tmp > .8
  mask <- LogicalNeuroVol(mask, drop_dim(spc))

  bv1 <- SparseNeuroVec(dat, spc, mask)
  bv2 <- bv1 + bv1
  bv3 <- bv1 * bv2
  bv4 <- bv3 - bv1
  expect_true(TRUE)
})

test_that("can add/subtract NeuroVol from NeuroVec", {
  dat <- array(rnorm(12*12*12*4), c(12,12,12,4))
  spc <- NeuroSpace(c(12,12,12,4))
  tmp <- rnorm(12*12*12)
  mask <- tmp > .8
  mask <- LogicalNeuroVol(mask, drop_dim(spc))

  bv1 <- SparseNeuroVec(dat, spc, mask)
  vol <- NeuroVol(array(rnorm(12*12*12),c(12,12,12)), space(mask))
  bv2 <- bv1 + vol
  bv3 <- bv2 - vol
  bv4 <- bv1 * vol
  expect_equal(sum(bv3), sum(bv1))
})



test_that("can concatenate a SparseNeuroVec", {
	dat <- array(0, c(12,12,12,4))
	spc <- NeuroSpace(c(12,12,12,4))
	tmp <- rnorm(12*12*12)

	mask <- tmp > .8
	mask <- LogicalNeuroVol(mask, drop_dim(spc))

	bv1 <- SparseNeuroVec(dat, spc, mask)

	bv2 <- concat(bv1, bv1)
	expect_true(inherits(bv2, "NeuroVec"))
	expect_equal(dim(bv2), c(12,12,12,8))

	bv3 <- concat(bv1,bv2, bv1, bv2)
	expect_true(inherits(bv3, "NeuroVec"))
	expect_equal(dim(bv3), c(12,12,12,24))

	#do.call(concat, list(bv1))
	#expect_equal(bv4[1,1,1,1],0)
})

test_that("can extract nonzero coords of SparseNeuroVec", {
  dat <- array(0, c(12,12,12,4))
  spc <- NeuroSpace(c(12,12,12,4))
  tmp <- rnorm(12*12*12)
  mask <- tmp > .8
  mask <- LogicalNeuroVol(mask, drop_dim(spc))

  bv1 <- SparseNeuroVec(dat, spc, mask)
  cds <- coords(bv1)
  ind <- indices(bv1)
  expect_equal(nrow(cds), sum(mask))
  expect_equal(ind, which(mask>0))
})

test_that("can write and read back an image vector", {
  mask <- read_vol(gmask)
  maskvec <- concat(mask,mask)
  fname <- paste(tempfile(), ".nii", sep="")
  write_vec(maskvec, fname)

  vec2 <- read_vec(fname)
  expect_true(all(maskvec == vec2))
  expect_true(identical(space(maskvec), space(vec2)))

})

test_that("can write and read back an image vector as a bigneurovec", {
  mask <- read_vol(gmask)
  maskvec <- concat(mask,mask)
  fname <- paste(tempfile(), ".nii", sep="")
  write_vec(maskvec, fname)

  vec2 <- read_vec(fname, mask=mask, mode="bigvec")

  v1 <- series(maskvec, which(mask > 0))
  v2 <- series(vec2, which(mask > 0))
  expect_true(identical(space(maskvec), space(vec2)))
  expect_true(all(v1==v2))

})

test_that("can write and read back an image vector as a filebacked neurovec", {
  mask <- read_vol(gmask)
  maskvec <- concat(mask,mask)
  fname <- paste(tempfile(), ".nii", sep="")
  write_vec(maskvec, fname)

  vec2 <- read_vec(fname, mask=mask, mode="filebacked")

  v1 <- series(maskvec, which(mask > 0))
  v2 <- series(vec2, which(mask > 0))
  expect_true(identical(space(maskvec), space(vec2)))
  expect_true(all(v1==v2))

})

test_that("can extract ROI from NeuroVec", {
  mask <- read_vol(gmask)
  maskvec <- concat(mask,mask,mask,mask)
  roi <- spherical_roi(mask, centroid=c(10,10,10), radius=8)
  cds <- coords(roi)
  roi <- series_roi(maskvec, cds)
  vec <- as(roi, "SparseNeuroVec")
  expect_equivalent(class(vec), "SparseNeuroVec")
  expect_equivalent(nrow(roi@coords), sum(vec@mask))

})

test_that("can extract ROI from SparseNeuroVec", {
  mask <- read_vol(gmask)
  maskvec <- concat(mask,mask,mask,mask)
  maskvec <- SparseNeuroVec(maskvec@.Data, space(maskvec), mask=mask)
  roi <- spherical_roi(mask, centroid=c(17,17,10), radius=8)
  cds <- coords(roi)
  roi <- series_roi(maskvec, cds)
  vec <- as(roi, "SparseNeuroVec")
  expect_equivalent(class(vec), "SparseNeuroVec")
  expect_equivalent(nrow(roi@coords), sum(vec@mask))

})

library(testthat)

test_that("NeuroVec constructor works correctly", {
  # Create a simple 4D array
  data <- array(1:24, dim = c(2, 2, 3, 2))

  # Create a NeuroSpace object
  space <- NeuroSpace(c(2, 2, 3, 2), c(1,1,1),
                      c(0, 0, 0))

  # Test NeuroVec constructor with a 4D array
  nv <- NeuroVec(data, space)
  expect_s4_class(nv, "DenseNeuroVec")
  expect_equal(nv@.Data, data)

  space2 <- NeuroSpace(c(2, 2, 3), c(1,1,1),
                      c(0, 0, 0))
  # Test NeuroVec constructor with a list of NeuroVols
  nv_list <- list(NeuroVol(array(1:12, c(2, 2, 3)), space2),
                  NeuroVol(array(13:24, c(2, 2, 3)), space2))
  nv2 <- NeuroVec(nv_list)
  expect_s4_class(nv2, "DenseNeuroVec")
  expect_equivalent(nv2@.Data, c(data))
})



# test.NeuroVec.roundtrip.io <- function() {
#   bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4)
#   template <- takeVolume(bvec,1)
#
#   mask.idx <- sort(sample(1:length(template), 1000))
#   vals <- rnorm(length(mask.idx))
#   bv <- NeuroVol(vals, space(template), indices=mask.idx)
#   fname <- paste(tempfile(), ".nii", sep="")
#   writeVolume(bv,fname)
#   bv2 <- loadVolume(fname)
#
#   expect_equal(dim(bv2), dim(bv))
#   expect_equalNumeric(trans(bv2), trans(bv), tol=.0001)
#
# }


