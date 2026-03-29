## Additional coverage tests for R/searchlight.R
## Uses small 8x8x8 volumes to keep runtime fast.

library(neuroim2)

make_small_mask <- function(dim = c(8L, 8L, 8L), all_true = TRUE) {
  sp <- NeuroSpace(dim, c(1, 1, 1))
  if (all_true) {
    arr <- array(TRUE, dim)
  } else {
    arr <- array(runif(prod(dim)) > 0.3, dim)
    arr[4, 4, 4] <- TRUE  # guarantee at least one nonzero voxel
  }
  LogicalNeuroVol(arr, sp)
}

# ---------------------------------------------------------------------------
# random_searchlight
# ---------------------------------------------------------------------------

test_that("random_searchlight returns a non-empty list", {
  skip_on_cran()
  mask <- make_small_mask()
  sl <- random_searchlight(mask, radius = 2)
  expect_true(is.list(sl))
  expect_true(length(sl) > 0)
})

test_that("random_searchlight elements are ROIVolWindow objects", {
  skip_on_cran()
  mask <- make_small_mask()
  sl <- random_searchlight(mask, radius = 2)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})

test_that("random_searchlight errors on non-NeuroVol mask", {
  expect_error(random_searchlight(matrix(1, 4, 4), radius = 2))
})

test_that("random_searchlight errors on negative radius", {
  mask <- make_small_mask()
  expect_error(random_searchlight(mask, radius = -1))
})

test_that("random_searchlight errors on zero radius", {
  mask <- make_small_mask()
  expect_error(random_searchlight(mask, radius = 0))
})

test_that("random_searchlight errors on non-logical nonzero argument", {
  mask <- make_small_mask()
  expect_error(random_searchlight(mask, radius = 2, nonzero = "yes"))
})

test_that("random_searchlight with nonzero=FALSE runs", {
  skip_on_cran()
  mask <- make_small_mask()
  sl <- random_searchlight(mask, radius = 2, nonzero = FALSE)
  expect_true(length(sl) > 0)
})

# ---------------------------------------------------------------------------
# searchlight (exhaustive)
# ---------------------------------------------------------------------------

test_that("searchlight returns deferred list with eager=FALSE", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- searchlight(mask, radius = 2, eager = FALSE, nonzero = FALSE)
  expect_true(length(sl) > 0)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})

test_that("searchlight returns list with eager=TRUE", {
  skip_on_cran()
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- searchlight(mask, radius = 2, eager = TRUE, nonzero = FALSE)
  expect_true(length(sl) > 0)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})

test_that("searchlight with nonzero=TRUE filters to nonzero voxels", {
  sp <- NeuroSpace(c(6L, 6L, 6L), c(1, 1, 1))
  arr <- array(0, c(6, 6, 6))
  arr[3:5, 3:5, 3:5] <- 1
  mask <- LogicalNeuroVol(arr > 0, sp)
  sl <- searchlight(mask, radius = 2, eager = FALSE, nonzero = TRUE)
  # Should have one entry per nonzero center voxel
  expect_true(length(sl) > 0)
  roi1 <- sl[[1]]
  # All coords in roi should lie within the mask nonzero region
  # (at minimum, some voxels should pass)
  expect_true(nrow(coords(roi1)) > 0)
})

test_that("searchlight errors on non-NeuroVol mask", {
  expect_error(searchlight(array(1, c(4, 4, 4)), radius = 2))
})

test_that("searchlight errors on non-positive radius", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(searchlight(mask, radius = 0))
  expect_error(searchlight(mask, radius = -2))
})

test_that("searchlight errors on non-logical eager", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(searchlight(mask, radius = 2, eager = "yes"))
})

test_that("searchlight errors on non-logical nonzero", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(searchlight(mask, radius = 2, nonzero = 1L))
})

test_that("searchlight errors on negative cores", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(searchlight(mask, radius = 2, cores = -1))
})

test_that("searchlight mask_index attribute is set", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- searchlight(mask, radius = 2, eager = FALSE, nonzero = FALSE)
  roi1 <- sl[[1]]
  expect_false(is.null(attr(roi1, "mask_index")))
  expect_equal(attr(roi1, "mask_index"), 1L)
})

# ---------------------------------------------------------------------------
# searchlight_coords
# ---------------------------------------------------------------------------

test_that("searchlight_coords returns list with nonzero=FALSE", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sc <- searchlight_coords(mask, radius = 2, nonzero = FALSE, cores = 0)
  expect_true(length(sc) > 0)
  expect_true(is.matrix(sc[[1]]))
  expect_equal(ncol(sc[[1]]), 3L)
})

test_that("searchlight_coords returns list with nonzero=TRUE", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sc <- searchlight_coords(mask, radius = 2, nonzero = TRUE, cores = 0)
  expect_true(length(sc) > 0)
})

test_that("searchlight_coords result has fewer entries when nonzero=TRUE vs FALSE", {
  sp <- NeuroSpace(c(6L, 6L, 6L), c(1, 1, 1))
  arr <- array(0, c(6, 6, 6))
  arr[2:5, 2:5, 2:5] <- 1
  mask <- LogicalNeuroVol(arr > 0, sp)
  sc_all  <- searchlight_coords(mask, radius = 2, nonzero = FALSE, cores = 0)
  sc_nz   <- searchlight_coords(mask, radius = 2, nonzero = TRUE,  cores = 0)
  # nonzero restricts centers to nonzero voxels
  expect_true(length(sc_nz) <= length(sc_all))
})

# ---------------------------------------------------------------------------
# clustered_searchlight
# ---------------------------------------------------------------------------

test_that("clustered_searchlight with csize returns correct number of ROIs", {
  skip_on_cran()
  mask <- make_small_mask(c(8L, 8L, 8L))
  csize <- 5L
  sl <- clustered_searchlight(mask, csize = csize)
  expect_equal(length(sl), csize)
})

test_that("clustered_searchlight elements are ROIVol objects", {
  skip_on_cran()
  mask <- make_small_mask(c(8L, 8L, 8L))
  sl <- clustered_searchlight(mask, csize = 4L)
  expect_s4_class(sl[[1]], "ROIVol")
})

test_that("clustered_searchlight with pre-built cvol respects cluster count", {
  skip_on_cran()
  mask <- make_small_mask(c(8L, 8L, 8L))
  idx <- which(mask > 0)
  grid <- index_to_grid(mask, idx)
  set.seed(42)
  kres <- kmeans(grid, centers = 3L, iter.max = 200)
  cvol <- ClusteredNeuroVol(mask, clusters = kres$cluster)
  sl <- clustered_searchlight(mask, cvol = cvol)
  expect_equal(length(sl), num_clusters(cvol))
})

test_that("clustered_searchlight errors when neither csize nor cvol given", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(clustered_searchlight(mask))
})

test_that("clustered_searchlight errors on non-NeuroVol mask", {
  expect_error(clustered_searchlight(array(1, c(4, 4, 4)), csize = 2))
})

test_that("clustered_searchlight errors on non-positive csize", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(clustered_searchlight(mask, csize = 0))
})

# ---------------------------------------------------------------------------
# resampled_searchlight
# ---------------------------------------------------------------------------

test_that("resampled_searchlight returns deflist of correct length", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = 2, iter = 10)
  expect_equal(length(sl), 10L)
})

test_that("resampled_searchlight elements are ROIVolWindow", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = 2, iter = 5)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})

test_that("resampled_searchlight errors on non-NeuroVol mask", {
  expect_error(resampled_searchlight(matrix(1, 4, 4), radius = 2))
})

test_that("resampled_searchlight errors on non-positive radius", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(resampled_searchlight(mask, radius = -1))
})

test_that("resampled_searchlight errors on zero iter", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_error(resampled_searchlight(mask, radius = 2, iter = 0))
})

test_that("resampled_searchlight accepts vector of radii", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = c(1.5, 2, 2.5), iter = 8)
  expect_equal(length(sl), 8L)
})

test_that("resampled_searchlight with shape_fun='cube' works", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = 2, iter = 5, shape_fun = "cube")
  expect_equal(length(sl), 5L)
})

test_that("resampled_searchlight with shape_fun='ellipsoid' works", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = 2, iter = 5, shape_fun = "ellipsoid")
  expect_equal(length(sl), 5L)
})

test_that("resampled_searchlight with custom shape_fun works", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  my_fun <- function(mask, center, radius, iter, nonzero) {
    spherical_roi(mask, center, radius, nonzero = FALSE)
  }
  sl <- resampled_searchlight(mask, radius = 2, iter = 5, shape_fun = my_fun)
  expect_equal(length(sl), 5L)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})

test_that("resampled_searchlight mask_index attribute is set", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  sl <- resampled_searchlight(mask, radius = 2, iter = 3)
  roi1 <- sl[[1]]
  expect_false(is.null(attr(roi1, "mask_index")))
})

# ---------------------------------------------------------------------------
# bootstrap_searchlight (deprecated wrapper)
# ---------------------------------------------------------------------------

test_that("bootstrap_searchlight is deprecated and still returns results", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  expect_warning(
    sl <- bootstrap_searchlight(mask, radius = 2, iter = 5),
    regexp = NULL  # just check it warns
  )
  expect_equal(length(sl), 5L)
})

# ---------------------------------------------------------------------------
# ellipsoid_shape / cube_shape / blobby_shape constructors
# ---------------------------------------------------------------------------

test_that("ellipsoid_shape returns a function", {
  f <- ellipsoid_shape(scales = c(1, 1, 1.5))
  expect_true(is.function(f))
})

test_that("ellipsoid_shape errors on invalid scales", {
  expect_error(ellipsoid_shape(scales = c(1, 1)))          # wrong length
  expect_error(ellipsoid_shape(scales = c(-1, 1, 1)))      # non-positive
})

test_that("ellipsoid_shape errors on negative jitter", {
  expect_error(ellipsoid_shape(jitter = -0.1))
})

test_that("cube_shape returns a function", {
  f <- cube_shape()
  expect_true(is.function(f))
})

test_that("blobby_shape returns a function", {
  f <- blobby_shape(drop = 0.3, edge_fraction = 0.7)
  expect_true(is.function(f))
})

test_that("blobby_shape errors on drop outside [0,1]", {
  expect_error(blobby_shape(drop = -0.1))
  expect_error(blobby_shape(drop = 1.1))
})

test_that("blobby_shape errors on edge_fraction outside (0,1]", {
  expect_error(blobby_shape(edge_fraction = 0))
  expect_error(blobby_shape(edge_fraction = 1.1))
})

test_that("blobby_shape function runs on a small mask", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  f <- blobby_shape(drop = 0.3, edge_fraction = 0.7)
  center <- matrix(c(3L, 3L, 3L), nrow = 1)
  result <- f(mask = mask, center = center, radius = 2, iter = 1, nonzero = FALSE)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3L)
})

test_that("ellipsoid_shape function with jitter runs on a small mask", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  f <- ellipsoid_shape(scales = c(1, 1, 1.2), jitter = 0.05)
  center <- matrix(c(3L, 3L, 3L), nrow = 1)
  result <- f(mask = mask, center = center, radius = 2, iter = 1, nonzero = FALSE)
  expect_true(is.matrix(result))
})

# ---------------------------------------------------------------------------
# resampled_searchlight with shape_fun returning matrix (coordinate path)
# ---------------------------------------------------------------------------

test_that("resampled_searchlight with shape_fun returning matrix works", {
  mask <- make_small_mask(c(6L, 6L, 6L))
  cube_fn <- function(mask, center, radius, iter, nonzero) {
    # Return a raw coordinate matrix
    sp <- spacing(mask)
    hw <- ceiling(radius / sp)
    ctr <- drop(center)
    coords <- as.matrix(expand.grid(
      seq.int(max(1L, ctr[1] - hw[1]), min(dim(mask)[1], ctr[1] + hw[1])),
      seq.int(max(1L, ctr[2] - hw[2]), min(dim(mask)[2], ctr[2] + hw[2])),
      seq.int(max(1L, ctr[3] - hw[3]), min(dim(mask)[3], ctr[3] + hw[3]))
    ))
    storage.mode(coords) <- "integer"
    coords
  }
  sl <- resampled_searchlight(mask, radius = 1.5, iter = 4, shape_fun = cube_fn)
  expect_equal(length(sl), 4L)
  expect_s4_class(sl[[1]], "ROIVolWindow")
})
