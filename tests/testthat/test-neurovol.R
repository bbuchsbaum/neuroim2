
gmask <- system.file("extdata", "global_mask.nii", package="neuroim2")

test_that("can construct NeuroVol from 3D array", {
	dat <- array(0, c(64,64,64))
	spc <- NeuroSpace(c(64,64,64))
	bv <- DenseNeuroVol(dat, spc)
	expect_true(!is.null(bv))
	expect_equal(bv[1,1,1], 0)
	expect_equal(bv[64,64,64], 0)
	expect_equal(dim(bv), c(64,64,64))
})

test_that("can construct NeuroVol from 1D array with indices", {
	dat <- rnorm(100)
	spc <- NeuroSpace(c(64,64,64))

	indices = seq(1,20000, length.out=100)
	bv <- DenseNeuroVol(dat, spc, indices=indices)
	expect_true(!is.null(bv))
	expect_equal(dim(bv), c(64,64,64))
})

test_that("can concatenate multiple NeuroVols to get a NeuroVec", {
	dat <- array(0, c(64,64,64))
	spc <- NeuroSpace(c(64,64,64))
	bv1 <- DenseNeuroVol(dat, spc)
	bv2 <- DenseNeuroVol(dat, spc)

	bv3 <- concat(bv1, bv2)
	expect_true(inherits(bv3, "NeuroVec"))
	expect_equal(dim(bv3), c(64,64,64,2))

	bv4 <- concat(bv1,bv2, bv1, bv2)
	expect_true(inherits(bv4, "NeuroVec"))
	expect_equal(dim(bv4), c(64,64,64,4))
	expect_equal(bv4[1,1,1,1],0)
})

test_that("can read a NIFTI vol", {
	vol <- read_vol(gmask)
	expect_true(!is.null(vol))
	expect_equal(dim(vol), c(64,64,25))
	expect_error(read_vol(gmask, index=5))
})




test_that("can perform arithmetic on NeuroVols", {
	vol1 <- read_vol(gmask)
	vol2 <- vol1

	vol.plus <- vol1 + vol2
	vol.check <- DenseNeuroVol(vol1@.Data + vol2@.Data, space(vol1))
	expect_equal(vol.plus, vol.check)

	vol.minus <- vol1 - vol2
	vol.check <- DenseNeuroVol(vol1@.Data - vol2@.Data, space(vol1))
	expect_equal(vol.minus, vol.check)

	vol.mult <- vol1 * vol2
	vol.check <- DenseNeuroVol(vol1@.Data * vol2@.Data, space(vol1))
	expect_equal(vol.mult, vol.check)

	vol.div <- (vol1+10)/vol2
	vol.check <- DenseNeuroVol((vol1@.Data+10) / vol2@.Data, space(vol1))
	expect_equal(vol.div, vol.check)

})

test_that("index_to_grid on NeuroVol checks out", {
	vol1 <- read_vol(gmask)

	i <- 65
	expect_equal(index_to_grid(vol1, i), matrix(c(1,2,1), nrow=1))
	expect_equal(index_to_grid(vol1, 1), matrix(c(1,1,1), nrow=1))
	expect_error(index_to_grid(vol1, 64*64*25 +1))
	expect_error(index_to_grid(vol1, -1))
})

test_that("can convert NeuroVol to LogicalNeuroVol", {
	vol1 <- read_vol(gmask)
	vol2 <- as(vol1, "LogicalNeuroVol")
	vol3 <- as.logical(vol1)
	expect_true(!is.null(vol2))
})


test_that("can compute mean of each slice with 'slices'", {
	vol1 <- read_vol(gmask)
	slices <- slices(vol1)
	mean.slice1 <- lapply(slices, mean)
	expect_equal(length(mean.slice1), 25)
	expect_equal(unlist(mean.slice1), apply(vol1, 3, mean))
})

# test_that({
# 	vol1 <- loadVolume("data/global_mask.nii")
# 	fname <- paste(tempfile(), ".nii", sep="")
# 	writeVolume(vol1, fname)
#
#
# 	vol2 <- loadVolume(fname)
# 	expect_true(all(vol1 == vol2))
# 	expect_equal(vol2@source@metaInfo@dataType, vol1@source@metaInfo@dataType)
#
# 	expect_true(identical(space(vol1), space(vol2)))
#
# 	fname <- paste(tempfile(), ".nii", sep="")
# 	writeVolume(vol1, fname, dataType="DOUBLE")
# 	vol3 <- loadVolume(fname)
# 	expect_equal(vol3@source@metaInfo@dataType, "DOUBLE")
#
# }


