test_that("NeuroSpace constructor creates valid 3D space", {
  sp <- make_space(c(64, 64, 32), spacing = c(3, 3, 4))
  expect_equal(dim(sp), c(64L, 64L, 32L))
  expect_equal(spacing(sp), c(3, 3, 4))
  expect_equal(origin(sp), c(0, 0, 0))
  expect_equal(ndim(sp), 3L)
})

test_that("NeuroSpace constructor creates valid 2D space", {
  sp <- NeuroSpace(c(128L, 128L), spacing = c(1.5, 1.5))
  expect_equal(dim(sp), c(128L, 128L))
  expect_equal(ndim(sp), 2L)
  expect_equal(spacing(sp), c(1.5, 1.5))
})

test_that("NeuroSpace constructor creates valid 4D space", {
  sp <- NeuroSpace(c(64L, 64L, 32L, 100L), spacing = c(3, 3, 4))
  expect_equal(dim(sp), c(64L, 64L, 32L, 100L))
  expect_equal(ndim(sp), 4L)
})

test_that("NeuroSpace from explicit affine trans", {
  aff <- diag(c(2, 2, 2, 1))
  aff[1:3, 4] <- c(-90, -126, -72)
  sp <- NeuroSpace(c(91, 109, 91), trans = aff)
  expect_equal(trans(sp), aff)
  # Actual voxel sizes are encoded in the affine, not the spacing slot
  expect_equal(voxel_sizes(trans(sp)), c(2, 2, 2))
})

test_that("NeuroSpace validates positive dimensions", {
  expect_error(NeuroSpace(c(0, 10, 10)))
  expect_error(NeuroSpace(c(-5, 10, 10)))
})

test_that("NeuroSpace validates positive spacing", {
  expect_error(NeuroSpace(c(10, 10, 10), spacing = c(-1, 1, 1)))
  expect_error(NeuroSpace(c(10, 10, 10), spacing = c(0, 1, 1)))
})

# ---- Coordinate transforms ----

test_that("index_to_grid and grid_to_index are inverses", {
  sp <- make_space(c(10, 10, 10))
  idx <- c(1, 50, 100, 500, 1000)
  grids <- index_to_grid(sp, idx)
  idx2 <- grid_to_index(sp, grids)
  expect_equal(idx, as.numeric(idx2))
})

test_that("grid_to_coord and coord_to_grid are inverses", {
  sp <- make_space(c(64, 64, 32), spacing = c(3, 3, 4), origin = c(-90, -126, -72))
  grid <- matrix(c(10, 20, 15), nrow = 1)
  coord <- grid_to_coord(sp, grid)
  grid2 <- coord_to_grid(sp, coord)
  expect_equal(as.numeric(grid), as.numeric(grid2), tolerance = 1e-6)
})

test_that("index_to_coord and coord_to_index are inverses", {
  sp <- make_space(c(20, 20, 20), spacing = c(2, 2, 2), origin = c(-20, -20, -20))
  idx <- c(1, 100, 500, 4000, 8000)
  coords <- index_to_coord(sp, idx)
  idx2 <- coord_to_index(sp, coords)
  expect_equal(idx, as.numeric(idx2))
})

test_that("coord transforms work with batch inputs", {
  sp <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  grids <- matrix(c(1,1,1, 5,5,5, 10,10,10), ncol = 3, byrow = TRUE)
  coords <- grid_to_coord(sp, grids)
  expect_equal(nrow(coords), 3)
  grids2 <- coord_to_grid(sp, coords)
  expect_equal(as.numeric(grids), as.numeric(grids2), tolerance = 1e-6)
})

# ---- Affine operations ----

test_that("trans() returns the transformation matrix", {
  sp <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  tx <- trans(sp)
  expect_true(is.matrix(tx))
  expect_equal(dim(tx), c(4, 4))
  # Diagonal should reflect spacing
  expect_equal(tx[1,1], 2)
  expect_equal(tx[2,2], 2)
  expect_equal(tx[3,3], 2)
})

test_that("inverse_trans() is the inverse of trans()", {
  sp <- make_space(c(10, 10, 10), spacing = c(2, 3, 4), origin = c(-10, -15, -20))
  tx <- trans(sp)
  itx <- inverse_trans(sp)
  product <- tx %*% itx
  expect_equal(product, diag(4), tolerance = 1e-6)
})

test_that("oblique affine round-trips through coord transforms", {
  aff <- matrix(c(1.5, 0.1, 0, -80,
                   0.1, 2.0, 0, -120,
                   0, 0, 2.5, -60,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(64, 64, 32), trans = aff)

  idx <- c(1, 100, 1000, 50000)
  coords <- index_to_coord(sp, idx)
  idx2 <- coord_to_index(sp, coords)
  expect_equal(idx, as.numeric(idx2), tolerance = 0.5)
})

# ---- Spatial queries ----

test_that("spacing() returns correct values", {
  sp <- make_space(c(10, 10, 10), spacing = c(1.5, 2.0, 2.5))
  expect_equal(spacing(sp), c(1.5, 2.0, 2.5))
})

test_that("origin() returns correct values", {
  sp <- make_space(c(10, 10, 10), origin = c(-90, -126, -72))
  expect_equal(origin(sp), c(-90, -126, -72))
})

test_that("bounds() returns min/max world coordinates", {
  sp <- make_space(c(10, 10, 10), spacing = c(1, 1, 1))
  b <- bounds(sp)
  expect_equal(nrow(b), 3)
  expect_equal(ncol(b), 2)
  for (i in 1:3) {
    expect_true(b[i, 1] <= b[i, 2])
  }
})

test_that("axes() returns an AxisSet", {
  sp <- make_space(c(10, 10, 10))
  ax <- axes(sp)
  expect_s4_class(ax, "AxisSet3D")
})

# ---- Dimension manipulation ----

test_that("add_dim extends space dimensionality", {
  sp <- make_space(c(10, 10, 10))
  sp4 <- add_dim(sp, 20)
  expect_equal(ndim(sp4), 4L)
  expect_equal(dim(sp4), c(10L, 10L, 10L, 20L))
})

test_that("drop_dim reduces space dimensionality", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 20L))
  sp3 <- drop_dim(sp4)
  expect_equal(ndim(sp3), 3L)
  expect_equal(dim(sp3), c(10L, 10L, 10L))
})

test_that("drop_dim with explicit dimnum works", {
  sp4 <- NeuroSpace(c(10L, 10L, 10L, 20L))
  sp3 <- drop_dim(sp4, 4)
  expect_equal(ndim(sp3), 3L)
  expect_equal(dim(sp3), c(10L, 10L, 10L))
})

test_that("drop_dim preserves spatial origin for 4D->3D", {
  sp3_orig <- make_space(c(10, 10, 10), spacing = c(2, 2, 2), origin = c(-10, -20, -30))
  sp4 <- add_dim(sp3_orig, 20)
  sp3 <- drop_dim(sp4)
  expect_equal(origin(sp3), c(-10, -20, -30))
  expect_equal(spacing(sp3), c(2, 2, 2))
})

test_that("add_dim then drop_dim round-trips for >3D", {
  sp <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  sp4 <- add_dim(sp, 50)
  sp3 <- drop_dim(sp4)
  expect_equal(dim(sp3), dim(sp))
  expect_equal(spacing(sp3), spacing(sp))
})

# ---- Equality ----

test_that("identical NeuroSpaces compare equal", {
  sp1 <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  sp2 <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  expect_equal(trans(sp1), trans(sp2))
  expect_equal(dim(sp1), dim(sp2))
})

test_that("different NeuroSpaces differ", {
  sp1 <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  sp2 <- make_space(c(10, 10, 10), spacing = c(3, 3, 3))
  expect_false(identical(trans(sp1), trans(sp2)))
})

# ---- Edge cases ----

test_that("single-voxel space has correct dimensions", {
  sp <- make_space(c(1, 1, 1))
  expect_equal(prod(dim(sp)), 1)
  expect_equal(ndim(sp), 3L)
  # grid_to_coord for grid (1,1,1): (1-1)*spacing + origin = (0,0,0)
  coord <- grid_to_coord(sp, matrix(c(1, 1, 1), nrow = 1))
  expect_equal(as.numeric(coord), c(0, 0, 0))
})

test_that("non-isotropic spacing is handled correctly", {
  sp <- make_space(c(20, 20, 10), spacing = c(1, 1, 3))
  expect_equal(spacing(sp), c(1, 1, 3))
  # grid_to_coord does (grid - 1) * spacing + origin
  # grid (10,10,5) -> (9, 9, 4) * (1, 1, 3) = (9, 9, 12)
  coord <- grid_to_coord(sp, matrix(c(10, 10, 5), nrow = 1))
  expect_equal(as.numeric(coord), c(9, 9, 12), tolerance = 1e-6)
})

test_that("centroid computes the mean world coordinate", {
  sp <- make_space(c(10, 10, 10), spacing = c(1, 1, 1))
  ctr <- centroid(sp)
  expect_equal(length(ctr), 3)
  # For uniform spacing from origin 0, centroid should be near the center
  expect_true(all(ctr > 0))
})

# ---- grid_to_grid ----

test_that("grid_to_grid permutes grid coordinates", {
  sp <- make_space(c(10, 10, 10), spacing = c(2, 2, 2))
  g <- matrix(c(1,2,3, 5,5,5), ncol = 3, byrow = TRUE)
  g2 <- grid_to_grid(sp, g)
  expect_equal(nrow(g2), 2)
  expect_equal(ncol(g2), 3)
})
