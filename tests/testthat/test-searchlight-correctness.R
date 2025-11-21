context("searchlight correctness and consistency")

test_that("searchlight eager == lazy and center_index matches voxel center", {
  sp <- NeuroSpace(c(5,5,5), spacing = c(1,1,1))
  mask <- LogicalNeuroVol(array(TRUE, dim = c(5,5,5)), sp)
  grid <- index_to_grid(mask, which(mask != 0))

  sl_lazy  <- searchlight(mask, radius = 1, eager = FALSE, nonzero = FALSE)
  sl_eager <- searchlight(mask, radius = 1, eager = TRUE,  nonzero = FALSE)

  for (i in 1:3) {
    rl <- sl_lazy[[i]]
    re <- sl_eager[[i]]
    expect_equal(coords(rl), coords(re))
    expect_equal(coords(rl)[rl@center_index, , drop = FALSE],
                 grid[i, , drop = FALSE])
    expect_equal(coords(re)[re@center_index, , drop = FALSE],
                 grid[i, , drop = FALSE])
  }
})

test_that("nonzero filtering shrinks searchlight neighborhoods", {
  sp <- NeuroSpace(c(5,5,5), spacing = c(1,1,1))
  arr <- array(FALSE, dim = c(5,5,5))
  arr[3,3,3] <- TRUE
  mask_sparse <- LogicalNeuroVol(arr, sp)

  sl_full <- searchlight(mask_sparse, radius = 2, eager = TRUE, nonzero = FALSE)
  sl_nz   <- searchlight(mask_sparse, radius = 2, eager = TRUE, nonzero = TRUE)

  expect_gt(nrow(coords(sl_full[[1]])), nrow(coords(sl_nz[[1]])))
  expect_equal(nrow(coords(sl_nz[[1]])), 1)
})

test_that("random_searchlight respects nonzero and tracks center_row_index", {
  sp <- NeuroSpace(c(5,5,5), spacing = c(1,1,1))
  arr <- array(FALSE, dim = c(5,5,5))
  arr[3,3,3] <- TRUE
  mask_sparse <- LogicalNeuroVol(arr, sp)

  set.seed(12)
  rs_full <- random_searchlight(mask_sparse, radius = 2, nonzero = FALSE)
  set.seed(12)
  rs_nz   <- random_searchlight(mask_sparse, radius = 2, nonzero = TRUE)

  expect_gt(nrow(coords(rs_full[[1]])), nrow(coords(rs_nz[[1]])))
  expect_equal(nrow(coords(rs_nz[[1]])), 1)
  expect_identical(attr(rs_nz[[1]], "center_row_index"), 1L)
})

test_that("resampled_searchlight honors nonzero and shape keywords", {
  sp <- NeuroSpace(c(5,5,5), spacing = c(1,1,1))
  arr <- array(FALSE, dim = c(5,5,5))
  arr[3,3,3] <- TRUE
  mask_sparse <- LogicalNeuroVol(arr, sp)

  set.seed(45)
  r_full <- resampled_searchlight(mask_sparse, radius = 2, iter = 2,
                                  nonzero = FALSE, shape_fun = "ellipsoid")
  set.seed(45)
  r_nz   <- resampled_searchlight(mask_sparse, radius = 2, iter = 2,
                                  nonzero = TRUE, shape_fun = "ellipsoid")

  expect_gt(nrow(coords(r_full[[1]])), nrow(coords(r_nz[[1]])))
  expect_equal(nrow(coords(r_nz[[1]])), 1)
  expect_equal(coords(r_nz[[1]])[r_nz[[1]]@center_index, , drop = FALSE],
               matrix(c(3,3,3), nrow = 1))
})

test_that("searchlight_coords applies nonzero filter", {
  sp <- NeuroSpace(c(5,5,5), spacing = c(1,1,1))
  arr <- array(FALSE, dim = c(5,5,5))
  arr[3,3,3] <- TRUE
  mask_sparse <- LogicalNeuroVol(arr, sp)

  coords_full <- searchlight_coords(mask_sparse, radius = 1, nonzero = FALSE)[[1]]
  coords_nz   <- searchlight_coords(mask_sparse, radius = 1, nonzero = TRUE)[[1]]

  expect_gt(nrow(coords_full), nrow(coords_nz))
  expect_equal(nrow(coords_nz), 1)
})
