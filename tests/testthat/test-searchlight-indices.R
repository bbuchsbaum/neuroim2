make_index_searchlight_mask <- function(
    dims = c(7L, 8L, 6L),
    voxel_spacing = c(1, 1, 1),
    active = NULL) {
  values <- array(FALSE, dims)
  if (is.null(active)) {
    values[2:6, 2:7, 2:5] <- TRUE
  } else {
    values[active] <- TRUE
  }
  LogicalNeuroVol(values, NeuroSpace(dims, spacing = voxel_spacing))
}

expect_indices_match_coords <- function(mask, radius, nonzero) {
  got <- searchlight_indices(mask, radius, nonzero = nonzero)
  coords <- searchlight_coords(mask, radius, nonzero = nonzero, cores = 0)

  expect_identical(length(got), length(coords))
  for (i in seq_along(got)) {
    converted <- index_to_grid(mask, got[[i]])
    storage.mode(converted) <- "integer"
    expect_identical(
      unname(converted),
      unname(coords[[i]]),
      info = paste("centre", i, "nonzero", nonzero)
    )
  }
  invisible(got)
}

test_that("searchlight_indices exposes an explicit index and domain contract", {
  mask <- make_index_searchlight_mask()
  got <- searchlight_indices(mask, radius = 2)

  expect_s3_class(got, "searchlight_indices")
  expect_true(is.list(got))
  expect_true(all(vapply(got, is.integer, logical(1))))
  expect_identical(attr(got, "center_indices"), which(mask != 0))
  expect_identical(attr(got, "space"), space(mask))
  expect_identical(attr(got, "radius"), 2)
  expect_identical(attr(got, "nonzero"), TRUE)
  expect_identical(length(got), sum(mask != 0))
})

test_that("searchlight_indices matches coordinate neighborhoods", {
  mask <- make_index_searchlight_mask()
  expect_indices_match_coords(mask, radius = 2, nonzero = FALSE)
  expect_indices_match_coords(mask, radius = 2, nonzero = TRUE)
})

test_that("searchlight_indices clips boundary neighborhoods", {
  mask <- make_index_searchlight_mask(
    dims = c(4L, 5L, 3L),
    active = c(1L, 2L, 5L, 6L, 21L)
  )
  got <- expect_indices_match_coords(mask, radius = 1, nonzero = FALSE)

  expect_identical(attr(got, "center_indices"), which(mask != 0))
  expect_true(all(unlist(got, use.names = FALSE) >= 1L))
  expect_true(all(unlist(got, use.names = FALSE) <= prod(dim(mask))))
})

test_that("searchlight_indices nonzero filtering excludes zero mask voxels", {
  mask <- make_index_searchlight_mask(
    dims = c(5L, 5L, 5L),
    active = c(63L, 64L)
  )
  unrestricted <- searchlight_indices(mask, radius = 1, nonzero = FALSE)
  restricted <- searchlight_indices(mask, radius = 1, nonzero = TRUE)
  allowed <- which(mask != 0)

  expect_true(any(!unlist(unrestricted, use.names = FALSE) %in% allowed))
  expect_true(all(unlist(restricted, use.names = FALSE) %in% allowed))
  expect_true(all(lengths(restricted) < lengths(unrestricted)))
})

test_that("searchlight_indices respects anisotropic millimetre spacing", {
  mask <- make_index_searchlight_mask(voxel_spacing = c(1, 2, 4))
  got <- expect_indices_match_coords(mask, radius = 4, nonzero = TRUE)
  centre <- which(mask != 0)[length(got) %/% 2L]
  centre_grid <- drop(index_to_grid(mask, centre))
  members <- index_to_grid(mask, got[[length(got) %/% 2L]])
  delta_mm <- sweep(sweep(members, 2, centre_grid), 2, c(1, 2, 4), "*")
  distances <- sqrt(rowSums(delta_mm^2))

  expect_true(all(distances <= 4))
})

test_that("searchlight_indices is deterministic and leaves future state alone", {
  mask <- make_index_searchlight_mask()
  before <- future::plan()
  first <- searchlight_indices(mask, radius = 2)
  second <- searchlight_indices(mask, radius = 2)
  after <- future::plan()

  expect_identical(first, second)
  expect_identical(after, before)
  expect_false("cores" %in% names(formals(searchlight_indices)))
})

test_that("searchlight_indices does not construct ROI objects", {
  mask <- make_index_searchlight_mask()
  local_mocked_bindings(
    spherical_roi = function(...) stop("unexpected ROI construction"),
    .new_roi_vol_window = function(...) stop("unexpected ROI construction"),
    .package = "neuroim2"
  )

  expect_no_error(got <- searchlight_indices(mask, radius = 2))
  expect_true(all(vapply(got, is.integer, logical(1))))
})

test_that("searchlight_indices preserves an empty spatial domain", {
  mask <- LogicalNeuroVol(
    array(FALSE, c(4L, 5L, 3L)),
    NeuroSpace(c(4L, 5L, 3L))
  )
  got <- searchlight_indices(mask, radius = 1)

  expect_s3_class(got, "searchlight_indices")
  expect_length(got, 0L)
  expect_identical(attr(got, "center_indices"), integer(0))
  expect_identical(attr(got, "space"), space(mask))
})

test_that("searchlight_indices validates its public arguments", {
  mask <- make_index_searchlight_mask()

  expect_error(searchlight_indices(as.array(mask), 2), "NeuroVol")
  expect_error(searchlight_indices(mask, 0), "positive finite")
  expect_error(searchlight_indices(mask, c(1, 2)), "single positive")
  expect_error(searchlight_indices(mask, Inf), "positive finite")
  expect_error(searchlight_indices(mask, 2, nonzero = NA), "TRUE or FALSE")
  expect_error(searchlight_indices(mask, 2, nonzero = 1), "TRUE or FALSE")

  coarse <- make_index_searchlight_mask(voxel_spacing = c(2, 2, 2))
  expect_error(searchlight_indices(coarse, 1), "radius' is too small")
})

test_that("searchlight_indices prints a compact summary", {
  mask <- make_index_searchlight_mask()
  got <- searchlight_indices(mask, radius = 2)

  expect_output(print(got), "<searchlight_indices>")
  expect_output(print(got), "Centres:")
  expect_output(print(got), "Radius: 2 mm")
})
