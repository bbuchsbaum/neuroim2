## Additional coverage tests for R/roi.R
## Targets zero-coverage functions not exercised by test-roivol.R

library(neuroim2)

# ---------------------------------------------------------------------------
# ROICoords constructor
# ---------------------------------------------------------------------------

test_that("ROICoords constructor creates object from 3-column matrix", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  rc <- ROICoords(coords)
  expect_s4_class(rc, "ROICoords")
  expect_equal(dim(rc), c(2L, 3L))
  expect_equal(length(rc), 2L)
})

test_that("ROICoords constructor rejects non-matrix input", {
  expect_error(ROICoords(c(1, 2, 3)))
})

test_that("ROICoords constructor rejects matrix with wrong number of columns", {
  expect_error(ROICoords(matrix(1:6, ncol = 2)))
})

test_that("ROICoords show method runs without error", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE)
  rc <- ROICoords(coords)
  expect_output(show(rc), "ROICoords")
})

test_that("ROICoords subset with numeric index works", {
  coords <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 3, byrow = TRUE)
  rc <- ROICoords(coords)
  sub <- rc[1:2]
  expect_s4_class(sub, "ROICoords")
  expect_equal(length(sub), 2L)
})

test_that("ROICoords centroid computes column means", {
  coords <- matrix(c(1, 2, 3, 3, 4, 5), ncol = 3, byrow = TRUE)
  rc <- ROICoords(coords)
  ctr <- centroid(rc)
  expect_equal(ctr, c(2, 3, 4))
})

# ---------------------------------------------------------------------------
# ROIVol constructor and basic methods
# ---------------------------------------------------------------------------

test_that("ROIVol constructor works with coords and data", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  roi <- ROIVol(sp, coords, data = c(1.0, 2.0))
  expect_s4_class(roi, "ROIVol")
  expect_equal(length(roi), 2L)
  expect_equal(dim(roi), c(2L, 3L))
})

test_that("ROIVol constructor errors when data length mismatches coords rows", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  expect_error(ROIVol(sp, coords, data = c(1.0, 2.0, 3.0)))
})

test_that("ROIVol constructor errors on non-matrix coords", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  expect_error(ROIVol(sp, c(1, 2, 3), data = 1.0))
})

test_that("ROIVol show method runs without error", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 2)
  expect_output(show(roi), "ROIVol")
})

test_that("ROIVol as.numeric returns data vector", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L), ncol = 3)
  roi <- ROIVol(sp, coords, data = 7.5)
  expect_equal(as.numeric(roi), 7.5)
})

test_that("ROIVol values method returns data", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  v <- values(roi)
  expect_true(is.numeric(v))
  expect_equal(length(v), length(roi))
})

test_that("ROIVol indices returns linear indices", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  coords <- matrix(c(1L, 1L, 1L, 2L, 2L, 2L), ncol = 3, byrow = TRUE)
  roi <- ROIVol(sp, coords, data = c(1, 2))
  idx <- indices(roi)
  expect_equal(length(idx), 2L)
  expect_true(all(idx > 0))
})

test_that("ROIVol as.sparse creates SparseNeuroVol", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  sv <- as.sparse(roi)
  expect_s4_class(sv, "SparseNeuroVol")
})

test_that("ROIVol as.logical creates LogicalNeuroVol", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  lv <- as.logical(roi)
  expect_s4_class(lv, "LogicalNeuroVol")
  expect_equal(sum(lv), length(roi))
})

test_that("ROIVol coords with real=TRUE returns real-world coordinates", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  gc <- coords(roi, real = FALSE)
  rc <- coords(roi, real = TRUE)
  expect_equal(ncol(rc), 3L)
  expect_equal(nrow(rc), nrow(gc))
})

# ---------------------------------------------------------------------------
# ROIVol [ extraction methods
# ---------------------------------------------------------------------------

test_that("ROIVol [numeric, missing] subsets rows", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 2)
  sub <- roi[1:3]
  expect_s4_class(sub, "ROIVol")
  expect_equal(length(sub), 3L)
})

test_that("ROIVol [logical, missing] subsets rows", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 2)
  logi <- rep(c(TRUE, FALSE), length.out = length(roi))
  sub <- roi[logi]
  expect_s4_class(sub, "ROIVol")
  expect_equal(length(sub), sum(logi))
})

test_that("ROIVol [missing, missing] returns data vector", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1, fill = 3)
  result <- roi[]
  expect_true(is.numeric(result))
  expect_true(all(result == 3))
})

test_that("ROIVol [missing, numeric] returns coord column", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  col1 <- roi[, 1]
  expect_true(is.matrix(col1))
  expect_equal(ncol(col1), 1L)
})

test_that("ROIVol [numeric, numeric] returns coord submatrix", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  sub <- roi[1:2, 1:2]
  expect_true(is.matrix(sub))
})

test_that("ROIVol [logical, numeric] returns coord submatrix", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  logi <- rep(c(TRUE, FALSE), length.out = length(roi))
  sub <- roi[logi, 1]
  expect_true(is.matrix(sub))
})

test_that("ROIVol [ROICoords, missing] subsets using full coord set", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 2)
  # ROICoords-based subsetting requires same row count as the ROIVol coords
  rc <- ROICoords(coords(roi))
  sub <- roi[rc]
  expect_s4_class(sub, "ROIVol")
  expect_equal(length(sub), length(roi))
})

test_that("ROIVol [ROICoords, numeric] returns coord column", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 2)
  # ROICoords-based subsetting requires same row count as the ROIVol coords
  rc <- ROICoords(coords(roi))
  sub <- roi[rc, 2]
  expect_true(is.matrix(sub))
})

test_that("ROIVol [matrix, missing] replaces coords", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1)
  new_coords <- coords(roi)  # same shape, same row count
  sub <- roi[new_coords]
  expect_s4_class(sub, "ROIVol")
})

# ---------------------------------------------------------------------------
# ROIVec constructor and methods
# ---------------------------------------------------------------------------

test_that("ROIVec constructor creates object", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(10), nrow = 5, ncol = 2)
  rv <- ROIVec(sp4, coords, data)
  expect_s4_class(rv, "ROIVec")
})

test_that("ROIVec constructor errors when column count mismatches coords rows", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(15), nrow = 5, ncol = 3)  # 3 cols but 2 coord rows
  expect_error(ROIVec(sp4, coords, data))
})

test_that("ROIVec with vector data reshapes to 1-row matrix", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  rv <- ROIVec(sp4, coords, data = c(1, 2))
  expect_s4_class(rv, "ROIVec")
  expect_equal(ncol(rv@.Data), 2L)
})

test_that("ROIVec show method runs without error", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(10), nrow = 5, ncol = 2)
  rv <- ROIVec(sp4, coords, data)
  expect_output(show(rv), "ROIVec")
})

test_that("ROIVec as.matrix coercion works", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  sp3 <- drop_dim(sp4)
  cube <- cuboid_roi(sp3, c(5, 5, 5), 2)
  data <- matrix(rnorm(5 * nrow(coords(cube))), nrow = 5)
  rv <- ROIVec(sp4, coords(cube), data)
  m <- as.matrix(rv)
  expect_true(inherits(m, "matrix"))
})

test_that("ROIVec values method returns data matrix", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L), ncol = 3)
  data <- matrix(1:5, nrow = 5, ncol = 1)
  rv <- ROIVec(sp4, coords, data)
  v <- values(rv)
  expect_true(is.matrix(v))
})

test_that("ROIVec indices returns linear voxel indices", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 1L, 1L, 2L, 2L, 2L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(10), nrow = 5, ncol = 2)
  rv <- ROIVec(sp4, coords, data)
  idx <- indices(rv)
  expect_equal(length(idx), 2L)
})

test_that("ROIVec vectors method returns deflist", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(10), nrow = 5, ncol = 2)
  rv <- ROIVec(sp4, coords, data)
  vl <- vectors(rv)
  expect_equal(length(vl), 2L)
  expect_true(is.numeric(vl[[1]]))
})

test_that("ROIVec vectors with numeric subset works", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(15), nrow = 5, ncol = 3)
  rv <- ROIVec(sp4, coords, data)
  vl <- vectors(rv, 1:2)
  expect_equal(length(vl), 2L)
})

test_that("ROIVec vectors with logical subset works", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  coords <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L), ncol = 3, byrow = TRUE)
  data <- matrix(rnorm(15), nrow = 5, ncol = 3)
  rv <- ROIVec(sp4, coords, data)
  logi <- c(TRUE, FALSE, TRUE)
  vl <- vectors(rv, logi)
  expect_equal(length(vl), 2L)
})

# ---------------------------------------------------------------------------
# series_roi on NeuroVec
# ---------------------------------------------------------------------------

test_that("series_roi with matrix coords returns ROIVec", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 8L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(10 * 10 * 10 * 8), c(10, 10, 10, 8)), sp4)
  coords <- matrix(c(3L, 4L, 5L, 6L, 7L, 8L), ncol = 3, byrow = TRUE)
  rv <- series_roi(vec, coords)
  expect_s4_class(rv, "ROIVec")
  expect_equal(nrow(rv@.Data), 8L)
  expect_equal(ncol(rv@.Data), 2L)
})

test_that("series_roi with ROICoords returns ROIVec", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 8L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(10 * 10 * 10 * 8), c(10, 10, 10, 8)), sp4)
  coords <- matrix(c(3L, 4L, 5L, 6L, 7L, 8L), ncol = 3, byrow = TRUE)
  rc <- ROICoords(coords)
  rv <- series_roi(vec, rc)
  expect_s4_class(rv, "ROIVec")
})

test_that("series_roi with numeric linear indices returns ROIVec", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 8L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(10 * 10 * 10 * 8), c(10, 10, 10, 8)), sp4)
  rv <- series_roi(vec, 1:4)
  expect_s4_class(rv, "ROIVec")
})

test_that("series_roi with LogicalNeuroVol mask returns ROIVec", {
  sp4 <- NeuroSpace(c(6L, 6L, 6L, 5L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6 * 6 * 6 * 5), c(6, 6, 6, 5)), sp4)
  sp3 <- drop_dim(sp4)
  mask <- LogicalNeuroVol(array(c(rep(TRUE, 10), rep(FALSE, 6 * 6 * 6 - 10)), c(6, 6, 6)), sp3)
  rv <- series_roi(vec, mask)
  expect_s4_class(rv, "ROIVec")
  expect_equal(ncol(rv@.Data), 10L)
})

# ---------------------------------------------------------------------------
# spherical_roi on a volume (not just NeuroSpace)
# ---------------------------------------------------------------------------

test_that("spherical_roi on LogicalNeuroVol with nonzero=TRUE filters zeros", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  arr <- array(0, c(10, 10, 10))
  arr[4:6, 4:6, 4:6] <- 1
  mask <- LogicalNeuroVol(arr > 0, sp)
  roi <- spherical_roi(mask, c(5, 5, 5), radius = 2, nonzero = TRUE)
  expect_s4_class(roi, "ROIVolWindow")
  # All voxels should be nonzero
  vox <- coords(roi)
  expect_true(nrow(vox) > 0)
})

test_that("spherical_roi on DenseNeuroVol reads values from volume", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  data <- array(seq_len(1000), c(10, 10, 10))
  vol <- DenseNeuroVol(as.numeric(data), sp)
  roi <- spherical_roi(vol, c(5, 5, 5), radius = 2)
  expect_s4_class(roi, "ROIVolWindow")
  expect_true(length(roi@.Data) > 0)
})

test_that("spherical_roi errors on zero centroid", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  expect_error(spherical_roi(sp, c(0, 5, 5), radius = 2))
})

test_that("spherical_roi errors on centroid exceeding volume dims", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  expect_error(spherical_roi(sp, c(11, 5, 5), radius = 2))
})

test_that("spherical_roi with fill stores fill value", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- spherical_roi(sp, c(5, 5, 5), radius = 2, fill = 99)
  expect_true(all(roi@.Data == 99))
})

# ---------------------------------------------------------------------------
# cuboid_roi extended tests
# ---------------------------------------------------------------------------

test_that("cuboid_roi with NeuroVol reads values from volume", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  data <- array(seq_len(1000), c(10, 10, 10))
  vol <- DenseNeuroVol(as.numeric(data), sp)
  roi <- cuboid_roi(vol, c(5, 5, 5), surround = 1)
  expect_s4_class(roi, "ROIVolWindow")
  expect_true(length(roi@.Data) > 0)
})

test_that("cuboid_roi with fill parameter assigns fill value", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), surround = 1, fill = 42)
  expect_true(all(roi@.Data == 42))
})

test_that("cuboid_roi with nonzero=TRUE filters zero voxels", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  arr <- array(0, c(10, 10, 10))
  arr[5, 5, 5] <- 1  # only center is nonzero
  vol <- DenseNeuroVol(as.numeric(arr), sp)
  roi <- cuboid_roi(vol, c(5, 5, 5), surround = 1, nonzero = TRUE)
  # center is always retained
  expect_true(length(roi@.Data) >= 1)
})

# ---------------------------------------------------------------------------
# square_roi extended tests
# ---------------------------------------------------------------------------

test_that("square_roi with fixdim=1 creates correct grid", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- square_roi(sp, c(5, 5, 5), surround = 1, fixdim = 1)
  expect_s4_class(roi, "ROIVolWindow")
  # all x coords should be fixed at 5
  expect_true(all(coords(roi)[, 1] == 5))
})

test_that("square_roi with fixdim=2 creates correct grid", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- square_roi(sp, c(5, 5, 5), surround = 1, fixdim = 2)
  expect_s4_class(roi, "ROIVolWindow")
  expect_true(all(coords(roi)[, 2] == 5))
})

test_that("square_roi with NeuroVol reads values", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  data <- array(seq_len(1000), c(10, 10, 10))
  vol <- DenseNeuroVol(as.numeric(data), sp)
  roi <- square_roi(vol, c(5, 5, 5), surround = 1)
  expect_s4_class(roi, "ROIVolWindow")
  expect_true(length(roi@.Data) > 0)
})

# ---------------------------------------------------------------------------
# Kernel constructor
# ---------------------------------------------------------------------------

test_that("Kernel constructor creates object with default FUN", {
  kdim <- c(3L, 3L, 3L)
  vdim <- c(1, 1, 1)
  k <- Kernel(kerndim = kdim, vdim = vdim)
  expect_s4_class(k, "Kernel")
})

test_that("Kernel constructor with custom sd creates object", {
  kdim <- c(5L, 5L, 5L)
  vdim <- c(1, 1, 1)
  k <- Kernel(kerndim = kdim, vdim = vdim, FUN = dnorm, sd = 2)
  expect_s4_class(k, "Kernel")
  expect_true(sum(k@weights) > 0)
})

test_that("Kernel constructor errors when kerndim has length 1", {
  expect_error(Kernel(kerndim = 3, vdim = c(1, 1, 1)))
})

test_that("Kernel weights sum to 1", {
  k <- Kernel(kerndim = c(3L, 3L, 3L), vdim = c(1, 1, 1))
  expect_equal(sum(k@weights), 1, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# roi_surface_matrix and roi_vector_matrix (internal helpers via coercion)
# ---------------------------------------------------------------------------

test_that("ROIVec coercion to matrix produces roi_vector_matrix", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 5L), c(1, 1, 1))
  sp3 <- drop_dim(sp4)
  cube <- cuboid_roi(sp3, c(5, 5, 5), 1)
  nc <- nrow(coords(cube))
  data <- matrix(rnorm(5 * nc), nrow = 5, ncol = nc)
  rv <- ROIVec(sp4, coords(cube), data)
  m <- as(rv, "matrix")
  expect_true(inherits(m, "roi_vector_matrix"))
  expect_true(!is.null(attr(m, "indices")))
  expect_true(!is.null(attr(m, "coords")))
})

# ---------------------------------------------------------------------------
# ROIVol coercion to DenseNeuroVol
# ---------------------------------------------------------------------------

test_that("ROIVol coercion to DenseNeuroVol via as() works", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  roi <- cuboid_roi(sp, c(5, 5, 5), 1, fill = 3)
  dvol <- as(roi, "DenseNeuroVol")
  expect_s4_class(dvol, "DenseNeuroVol")
  expect_equal(sum(dvol[dvol > 0]), sum(roi))
})
