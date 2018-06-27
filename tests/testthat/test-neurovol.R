
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

test_that("can construct NeuroVol from 1D vector", {

  spc <- NeuroSpace(c(64,64,64))
  dat <- rnorm(prod(dim(spc)))
  bv <- DenseNeuroVol(dat, spc)
  expect_true(!is.null(bv))
  expect_equal(dim(bv), c(64,64,64))
})


test_that("can construct NeuroVol from matrix with 1 row or column", {

  spc <- NeuroSpace(c(64,64,64))
  dat <- as.matrix(rnorm(prod(dim(spc))))
  bv <- DenseNeuroVol(dat, spc)
  bv2 <- DenseNeuroVol(t(dat), spc)
  expect_equal(bv, bv2)
})

test_that("DenseNeuroVol throws error when space has wrong dimensions", {

  spc <- NeuroSpace(c(64,64))
  dat <- as.matrix(rnorm(prod(64*64*64)))
  expect_error(DenseNeuroVol(dat, spc))

})

test_that("can convert NeuroVol to LogicalNeuroVol", {

  spc <- NeuroSpace(c(64,64,64))
  dat <- rnorm(prod(dim(spc)))
  tot <- sum(dat>0)
  bv <- DenseNeuroVol(dat, spc)
  bvlog <- as.logical(bv>0)
  expect_equal(tot, sum(bvlog))
})


test_that("NeuroVol with mismatching data and space dimensions throw error", {

  spc <- NeuroSpace(c(64,64,64))
  dat <- array(0, c(64,64,63))
  expect_error(DenseNeuroVol(dat, spc))

})

test_that("can subset a NeuroVol with an ROIVol", {

  spc <- NeuroSpace(c(64,64,64))
  dat <- array(0, c(64,64,64))
  vol <- DenseNeuroVol(dat, spc)
  roi <- spherical_roi(vol, c(32,32,32), 3)
  vals <- vol[roi]
  expect_equal(length(vals), nrow(coords(roi)))
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
	vol4 <- as.mask(vol1)
	vol5 <- as.mask(vol1, indices=which(vol1>0))
	expect_true(!is.null(vol2))
	expect_equal(sum(vol3), sum(vol4))
	expect_equal(sum(vol4), sum(vol5))
})

test_that("can convert NeuroVol to SparseNeuroVol", {
  vol1 <- read_vol(gmask)
  svol1 <- as.sparse(vol1, mask=which(vol1>0))
  svol2 <- as.sparse(vol1, mask=as.logical(vol1))
  expect_equal(svol1,svol2)

})

test_that("data and indices must have same length to contruct SparseNeuroVol", {
  vol1 <- read_vol(gmask)
  indices <- which(vol1>0)
  dat <- vol1[indices]
  dat <- c(0,dat)
  expect_error(SparseNeuroVol(data=dat, indices=indices, space=space(vol1)))
})

test_that("can convert SparseNeuroVol to an array or vector", {
  vol1 <- read_vol(gmask)
  svol1 <- as.sparse(vol1, mask=which(vol1>0))

  arr <- as.array(svol1)
  vec <- as.numeric(svol1)
  expect_equal(dim(arr), dim(vol1))
  expect_equal(sum(vec), sum(svol1))
})

test_that("can map a kernel over a NeuroVol", {
  vol1 <- read_vol(gmask)
  kern <- Kernel(c(3,3,3), vdim=spacing(vol1))
  vol2 <- mapf(vol1, kern)
  expect_equal(space(vol1), space(vol2))
})

test_that("can map values over a NeuroVol", {
  vol1 <- read_vol(gmask)
  lookup <- list("0"=2,"1"=5)
  x <- map_values(vol1, lookup)
  expect_equal(names(table(x)), c("2", "5"))
  x2 <- map_values(vol1, rbind(c(0,2), c(1,5)))
  expect_equal(names(table(x2)), c("2", "5"))

})




test_that("can compute mean of each slice with 'slices'", {
	vol1 <- read_vol(gmask)
	slices <- slices(vol1)
	mean.slice1 <- lapply(slices, mean)
	expect_equal(length(mean.slice1), 25)
	expect_equal(unlist(mean.slice1), apply(vol1, 3, mean))
})

test_that("can write and read back an image vol", {
  mask <- read_vol(gmask)
	fname <- paste(tempfile(), ".nii", sep="")
	write_vol(mask, fname)

	vol2 <- read_vol(fname)
	expect_true(all(mask == vol2))
	expect_true(identical(space(mask), space(vol2)))

})


test_that("can construct a NeuroSlice", {
  slice <- NeuroSlice(matrix(0,10,10), NeuroSpace(c(10,10)))
  expect_true(!is.null(slice))
})

test_that("can construct a NeuroSlice and extract index", {
  slice <- NeuroSlice(matrix(0,10,10), NeuroSpace(c(10,10)))
  pix <- grid_to_index(slice, c(1,1))
  expect_equal(pix, 1)
  pix <- index_to_grid(slice, 1)
  expect_equal(as.vector(pix), c(1,1))
})

test_that("NeuroSlice with wrong input space throws error", {
  expect_error(NeuroSlice(matrix(0,10,10), NeuroSpace(c(10,10,10))))

})

test_that("NeuroSlice with mismatching data dim throws error", {
  expect_error(NeuroSlice(array(0,c(10,10,10)), NeuroSpace(c(10,10))))
})

test_that("Can construct a ClusteredNeuroVol", {
  mask <- NeuroVol(array(1,c(10,10,10)), NeuroSpace(c(10,10,10)))
  mask.idx <- which(mask != 0)
	grid <- index_to_coord(mask, mask.idx)
  vox <- index_to_grid(mask, mask.idx)
  kres <- kmeans(grid, centers=50, iter.max=500)
  kvol <- ClusteredNeuroVol(mask, kres$cluster)
  expect_true(!is.null(kvol))
  expect_equal(length(unique(kvol@clusters)), 50)
  expect_equal(length(kvol@cluster_map), 50)
  expect_equal(num_clusters(kvol), 50)
})

test_that("Can partition a NeuroVol into a set of ROIVols", {
  mask <- NeuroVol(array(1,c(10,10,10)), NeuroSpace(c(10,10,10)))
  mask.idx <- which(mask != 0)
  grid <- index_to_coord(mask, mask.idx)
  vox <- index_to_grid(mask, mask.idx)
  kres <- kmeans(grid, centers=50, iter.max=500)
  kvol <- ClusteredNeuroVol(mask, kres$cluster)
  p <- split_clusters(mask, kvol)
  expect_equal(50, length(p))
})



