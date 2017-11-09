
gen_dat <- function(d1=64,d2=64,d3=64,d4=4) {
   dat <- array(0, c(d1,d2,d3,d4))
   spc <- NeuroSpace(c(d1,d2,d3,d4))
   DenseNeuroVec(dat, spc)
}

test_that("can construct a DenseNeuroVec", {
	bv <- gen_dat(64,64,64,4)
	expect_true(!is.null(bv))
	expect_equal(bv[1,1,1,1], 0)
	expect_equal(bv[64,64,64,4], 0)
	expect_equal(dim(bv), c(64,64,64,4))
})

test_that("can concatenate two or more NeuroVecs", {
	bv1 <- gen_dat()
	bv2 <- gen_dat()

	bv3 <- concat(bv1, bv2)
	expect_true(inherits(bv3, "NeuroVec"))
	expect_equal(dim(bv3), c(64,64,64,8))

	bv4 <- concat(bv1,bv2, bv1, bv3)
	expect_true(inherits(bv4, "NeuroVec"))
	expect_equal(dim(bv4), c(64,64,64,20))
	expect_equal(bv4[1,1,1,1],0)

})

test_that("can extract a single volume from a NeuroVec", {
	bv1 <- gen_dat()
	bv2 <- gen_dat()

	bv3 <- concat(bv1, bv2)

	vol1 <- bv3[[1]]
	expect_equal(dim(vol1), c(64,64,64))
})


test_that("can extract a sub vector from a NeuroVec", {
  bv1 <- gen_dat()
  vec1 <- sub_vector(bv1, 1:2)
  expect_true(inherits(vec1, "NeuroVec"))
  expect_equal(dim(vec1), c(64,64,64,2))
})

test_that("can map over each volume in a NeuroVec", {
	bv1 <- gen_dat()
	mean.vol1 <- lapply(vols(bv1), mean)
	expect_equal(length(mean.vol1), 4)
	expect_equal(unlist(mean.vol1), apply(bv1, 4, mean))
})


test_that("can extract multiple series from a NeuroVec", {
	bv1 <- gen_dat()
	expect_equal(series(bv1, 1,1,1), bv1[1,1,1,])
  mat <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))

	r1 <- apply(mat, 1, function(i) { series(bv1, i[1], i[2], i[3]) })
	r2 <- series(bv1, mat)
	expect_equal(r1, r2)
})


test_that("can convert NeuroVec to matrix", {
	bv1 <- gen_dat(5,5,5,5)
	mat <- as.matrix(bv1)

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
# 	mask.idx <- sort(sample(1:length(template), 1000))
# 	vals <- rnorm(length(mask.idx))
# 	bv <- BrainVolume(vals, space(template), indices=mask.idx)
# 	fname <- paste(tempfile(), ".nii", sep="")
# 	writeVolume(bv,fname)
# 	bv2 <- loadVolume(fname)
#
# 	expect_equal(dim(bv2), dim(bv))
# 	expect_equalNumeric(trans(bv2), trans(bv), tol=.0001)
#
# }

test_that("can convert dense NeuroVec to sparse", {
	bv1 <- gen_dat()
	svec <- as.sparse(bv1, c(1,100,1000))
	expect_equal(dim(svec), c(64,64,64,4))

})

test_that("can construct a SparseNeuroVol", {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- NeuroSpace(c(64,64,64,4))
	tmp <- rnorm(64*64*64)
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
  dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
  spc <- NeuroSpace(c(64,64,64,4))
  tmp <- rnorm(64*64*64)
  mask <- tmp > .8
  mask <- LogicalNeuroVol(mask, drop_dim(spc))

  bv1 <- SparseNeuroVec(dat, spc, mask)
  bv2 <- bv1 + bv1
  bv3 <- bv1 * bv2
  bv4 <- bv3 - bv1
})



test_that("can concatenate a SparseNeuroVec", {
	dat <- array(0, c(64,64,64,4))
	spc <- NeuroSpace(c(64,64,64,4))
	tmp <- rnorm(64*64*64)
	mask <- tmp > .8
	mask <- LogicalNeuroVol(mask, drop_dim(spc))

	bv1 <- SparseNeuroVec(dat, spc, mask)

	bv2 <- concat(bv1, bv1)
	expect_true(inherits(bv2, "NeuroVec"))
	expect_equal(dim(bv2), c(64,64,64,8))

	bv3 <- concat(bv1,bv2, bv1, bv2)
	expect_true(inherits(bv3, "NeuroVec"))
	expect_equal(dim(bv3), c(64,64,64,24))

	do.call(concat, list(bv1))
	#expect_equal(bv4[1,1,1,1],0)
})



# test.NeuroVec.roundtrip.io <- function() {
#   bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4)
#   template <- takeVolume(bvec,1)
#
#   mask.idx <- sort(sample(1:length(template), 1000))
#   vals <- rnorm(length(mask.idx))
#   bv <- BrainVolume(vals, space(template), indices=mask.idx)
#   fname <- paste(tempfile(), ".nii", sep="")
#   writeVolume(bv,fname)
#   bv2 <- loadVolume(fname)
#
#   expect_equal(dim(bv2), dim(bv))
#   expect_equalNumeric(trans(bv2), trans(bv), tol=.0001)
#
# }


