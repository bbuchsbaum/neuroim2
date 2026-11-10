# The six coordinate conversions are defined on NeuroSpace. These tests pin
# that they also work
#   * for a single point (#54),
#   * on a volume or a vector, not only on its space (#58),
#   * whether the index or coordinate is stored as integer or double,
#   * for a data frame of points.

coord_fixture <- function() {
  tx <- rbind(c(-3, 0, 0, 90),
              c(0, 3, 0, -126),
              c(0, 0, 4, -72),
              c(0, 0, 0, 1))
  sp <- NeuroSpace(c(6L, 7L, 5L), trans = tx)
  vol <- NeuroVol(array(seq_len(prod(dim(sp))), dim(sp)), sp)
  vec <- NeuroVec(array(0, c(6, 7, 5, 2)), add_dim(sp, 2))
  list(sp = sp, vol = vol, vec = vec,
       mask = as.mask(vol > 100),
       sparse = as.sparse(vol, which(vol > 100)))
}

test_that("index_to_coord() takes a single index (#54)", {
  f <- coord_fixture()
  one <- index_to_coord(f$sp, 5)
  expect_true(is.matrix(one))
  expect_equal(dim(one), c(1L, 3L))
  expect_equal(one, grid_to_coord(f$sp, index_to_grid(f$sp, 5)))
  # one index is the first row of two
  expect_equal(one, index_to_coord(f$sp, c(5, 6))[1, , drop = FALSE])
  # integer and double agree
  expect_equal(index_to_coord(f$sp, 5L), one)
  expect_equal(index_to_coord(f$sp, 5:6), index_to_coord(f$sp, c(5, 6)))
})

test_that("index_to_coord() works on volumes and vectors (#58)", {
  f <- coord_fixture()
  idx <- c(1, 17, 100, 210)
  want <- grid_to_coord(f$sp, arrayInd(idx, dim(f$sp)))
  for (nm in c("vol", "vec", "mask", "sparse")) {
    expect_equal(index_to_coord(f[[nm]], idx), want, info = nm)
    expect_equal(index_to_coord(f[[nm]], as.integer(idx)), want, info = nm)
    expect_equal(index_to_coord(f[[nm]], 100), want[3, , drop = FALSE], info = nm)
    expect_equal(index_to_coord(f[[nm]], 100L), want[3, , drop = FALSE], info = nm)
  }
})

test_that("the peak of a map converts to millimetres in one call", {
  f <- coord_fixture()
  peak <- which.max(f$vol)
  expect_equal(as.vector(index_to_coord(f$vol, peak)),
               as.vector(grid_to_coord(f$sp, dim(f$sp))))
})

test_that("every conversion accepts a volume or vector in place of its space", {
  f <- coord_fixture()
  idx <- c(3, 50, 199)
  grid <- index_to_grid(f$sp, idx)
  world <- grid_to_coord(f$sp, grid)

  for (nm in c("vol", "vec", "mask", "sparse")) {
    x <- f[[nm]]
    # a 4-D object has a 4-D grid: these indices fall in its first volume
    grid_x <- if (nm == "vec") cbind(grid, 1) else grid
    expect_equal(index_to_grid(x, idx), grid_x, info = nm)
    expect_equal(index_to_grid(x, as.integer(idx)), grid_x, info = nm)
    expect_equal(grid_to_coord(x, grid), world, info = nm)
    expect_equal(grid_to_index(x, grid), idx, info = nm)
    expect_equal(coord_to_grid(x, world), grid, tolerance = 1e-5, info = nm)
    expect_equal(coord_to_index(x, world), idx, info = nm)

    # a single point as a plain vector
    expect_equal(as.vector(grid_to_coord(x, grid[2, ])), world[2, ], info = nm)
    expect_equal(grid_to_index(x, grid[2, ]), idx[2], info = nm)
    expect_equal(as.vector(coord_to_grid(x, world[2, ])), grid[2, ], tolerance = 1e-5, info = nm)
    expect_equal(coord_to_index(x, world[2, ]), idx[2], info = nm)
  }
})

test_that("integer-typed grids and coordinates are accepted", {
  f <- coord_fixture()
  g <- rbind(c(2L, 3L, 4L), c(6L, 7L, 5L))
  expect_equal(grid_to_coord(f$vol, g), grid_to_coord(f$sp, g + 0))
  expect_equal(grid_to_index(f$vol, g), grid_to_index(f$sp, g + 0))
  expect_equal(grid_to_coord(f$vol, c(2L, 3L, 4L)), grid_to_coord(f$sp, c(2, 3, 4)))
  w <- c(90L, -126L, -72L)   # the centre of voxel (1, 1, 1)
  expect_equal(coord_to_index(f$vol, w), 1)
  expect_equal(as.vector(coord_to_grid(f$vol, w)), c(1, 1, 1), tolerance = 1e-5)
})

test_that("a data frame of points is one point per row", {
  f <- coord_fixture()
  grid <- index_to_grid(f$sp, c(3, 50, 199))
  world <- grid_to_coord(f$sp, grid)
  wdf <- data.frame(x = world[, 1], y = world[, 2], z = world[, 3])
  gdf <- data.frame(i = grid[, 1], j = grid[, 2], k = grid[, 3])

  for (x in list(f$sp, f$vol, f$vec)) {
    expect_equal(coord_to_index(x, wdf), c(3, 50, 199))
    expect_equal(coord_to_grid(x, wdf), grid, tolerance = 1e-5)
    expect_equal(grid_to_coord(x, gdf), world)
    expect_equal(grid_to_index(x, gdf), c(3, 50, 199))
  }

  expect_error(coord_to_index(f$vol, cbind(wdf, label = "a")), "numeric columns")
})

test_that("points with the wrong number of columns are named as such", {
  f <- coord_fixture()
  expect_error(coord_to_index(f$sp, matrix(0, 4, 2)), "3 columns")
  expect_error(coord_to_grid(f$vol, matrix(0, 4, 2)), "3 columns")
  expect_error(grid_to_coord(f$vol, matrix(1, 4, 4)), "3 columns")
  # three points given as columns: the commonest slip
  expect_error(coord_to_index(f$sp, matrix(0, 3, 5)), "t\\(coords\\)")
})
