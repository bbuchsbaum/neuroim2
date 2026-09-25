library(testthat)

marker_xy <- function(p, marker = 1) {
  hit <- p$data[p$data$value == marker, c("x", "y"), drop = FALSE]
  expect_equal(nrow(hit), 1L)
  unname(as.numeric(hit[1L, ]))
}


test_that("generic plot and plot_montage preserve axial x/y landmark placement", {
  dims <- c(5L, 7L, 3L)
  arr <- array(0, dims)
  arr[2L, 6L, 2L] <- 1

  identity_vol <- NeuroVol(arr, NeuroSpace(dims))
  generic <- plot(identity_vol, zlevels = 2L)
  montage <- plot_montage(identity_vol, zlevels = 2L, ncol = 1L, range = "data")

  # Voxel axis 1 is screen x and voxel axis 2 is screen y. The old generic
  # path transposed these to (6, 2), rotating anterior toward screen-right.
  expect_equal(marker_xy(generic), c(1, 5))
  expect_equal(marker_xy(montage), c(1, 5))

  # A LAS-like x flip changes only left/right placement; anterior remains up.
  las <- diag(c(-1, 1, 1, 1))
  las[1, 4] <- dims[[1L]] - 1L
  las_vol <- NeuroVol(arr, NeuroSpace(dims, trans = las))
  generic_las <- plot(las_vol, zlevels = 2L)
  montage_las <- plot_montage(las_vol, zlevels = 2L, ncol = 1L, range = "data")
  las_slice <- slice(las_vol, 2L, along = 3L)
  slice_las <- plot(las_slice)
  slice_list_las <- plot_montage(list(las_slice), ncol = 1L, range = "data")

  expect_equal(marker_xy(generic_las), c(3, 5))
  expect_equal(marker_xy(montage_las), c(3, 5))
  expect_equal(marker_xy(slice_las), c(3, 5))
  expect_equal(marker_xy(slice_list_las), c(3, 5))
})

test_that("plot_montage orients a permuted native plane anatomically", {
  dims <- c(5L, 6L, 7L)
  arr <- array(0, dims)
  arr[2L, 5L, 3L] <- 1

  # Native axes increase toward Anterior, Superior, and Right (ASR).
  affine <- diag(4)
  affine[1:3, 1:3] <- matrix(
    c(0, 1, 0,
      0, 0, 1,
      1, 0, 0),
    nrow = 3L
  )
  vol <- NeuroVol(arr, NeuroSpace(dims, trans = affine))

  # Holding native axis 3 fixed produces a sagittal plane. Its horizontal
  # direction is P->A (native axis 1) and vertical direction I->S (axis 2).
  p <- plot_montage(vol, zlevels = 3L, along = 3L, ncol = 1L, range = "data")
  expect_equal(marker_xy(p), c(1, 4))
  expect_warning(ggplot2::ggplot_build(p), NA)
})

test_that("plot_ortho finds anatomical planes after voxel-axis permutation", {
  dims <- c(5L, 6L, 7L)
  marker <- c(2L, 5L, 3L)
  arr <- array(0, dims)
  arr[matrix(marker, nrow = 1L)] <- 1
  affine <- diag(4)
  affine[1:3, 1:3] <- matrix(
    c(0, 1, 0,
      0, 0, 1,
      1, 0, 0),
    nrow = 3L
  )
  vol <- NeuroVol(arr, NeuroSpace(dims, trans = affine))

  panels <- plot_ortho(vol, coord = marker, draw = FALSE, assemble = FALSE)
  expect_named(panels, c("axial", "coronal", "sagittal"))
  expect_equal(marker_xy(panels$axial), c(2, 1))
  expect_equal(marker_xy(panels$coronal), c(2, 4))
  expect_equal(marker_xy(panels$sagittal), c(1, 4))

  sides <- c("left", "right", "top", "bottom")
  expect_equal(tile_orientation_letters(panels$axial)[sides],
               c(left = "L", right = "R", top = "A", bottom = "P"))
  expect_equal(tile_orientation_letters(panels$coronal)[sides],
               c(left = "L", right = "R", top = "S", bottom = "I"))
  expect_equal(tile_orientation_letters(panels$sagittal)[sides],
               c(left = "P", right = "A", top = "S", bottom = "I"))

  # Each view is labelled with the world coordinate of its plane. The affine
  # maps native (i, j, k) to world (k, i, j) - 1, so the marker sits at
  # world (2, 1, 4).
  expect_world_label(tile_slice_label(panels$axial), "z = 4")
  expect_world_label(tile_slice_label(panels$coronal), "y = 1")
  expect_world_label(tile_slice_label(panels$sagittal), "x = 2")

  # The selected landmark and both crosshair arms must meet at the same
  # display coordinate in every plane, with a gap that leaves the landmark
  # itself uncovered.
  for (panel in panels) {
    xy <- marker_xy(panel)
    segs <- tile_crosshair(panel)
    vertical <- segs[segs$x == segs$xend, , drop = FALSE]
    horizontal <- segs[segs$y == segs$yend, , drop = FALSE]
    expect_gt(nrow(vertical), 0L)
    expect_gt(nrow(horizontal), 0L)
    expect_true(all(vertical$x == xy[[1L]]))
    expect_true(all(horizontal$y == xy[[2L]]))
    expect_false(any(vertical$y <= xy[[2L]] & vertical$yend >= xy[[2L]]))
    expect_false(any(horizontal$x <= xy[[1L]] & horizontal$xend >= xy[[1L]]))
  }
})

test_that("oblique plotting keeps native pixels on a regular warning-free grid", {
  dims <- c(8L, 10L, 6L)
  arr <- array(0, dims)
  arr[2L, 8L, 3L] <- 1
  affine <- diag(4)
  affine[1, 2] <- 0.25
  affine[2, 1] <- -0.15
  vol <- NeuroVol(arr, NeuroSpace(dims, trans = affine))

  generic <- plot(vol, zlevels = 3L)
  montage <- plot_montage(vol, zlevels = 3L, ncol = 1L, range = "data")
  ortho <- plot_ortho(vol, coord = c(2L, 8L, 3L), draw = FALSE, assemble = FALSE)

  expect_warning(ggplot2::ggplot_build(generic), NA)
  expect_warning(ggplot2::ggplot_build(montage), NA)
  expect_length(ortho, 3L)
  for (panel in ortho) {
    expect_warning(ggplot2::ggplot_build(panel), NA)
  }

  # The affine is oblique, but a native slice is still a regular voxel lattice.
  # Display coordinates therefore have one interval per native in-plane axis.
  expect_length(unique(round(diff(sort(unique(montage$data$x))), 12)), 1L)
  expect_length(unique(round(diff(sort(unique(montage$data$y))), 12)), 1L)
})

test_that("overlay and registration QC plots share the orientation transform", {
  dims <- c(6L, 8L, 4L)
  affine <- diag(4)
  affine[1, 2] <- 0.2
  affine[2, 1] <- -0.1
  sp <- NeuroSpace(dims, trans = affine)
  arr <- array(seq_len(prod(dims)), dims)
  bg <- NeuroVol(arr, sp)
  ov <- NeuroVol(rev(arr), sp)
  edges <- NeuroVol(array(as.numeric(arr %% 7L == 0L), dims), sp)

  overlay <- plot_overlay(
    bg, ov, zlevels = 2L, draw = FALSE, assemble = FALSE
  )[[1L]]
  checker <- plot_checkerboard(bg, ov, zlevels = 2L, draw = FALSE, assemble = FALSE)[[1L]]
  edge <- plot_edge_overlay(bg, edges, edges, zlevels = 2L, draw = FALSE,
                            assemble = FALSE)[[1L]]

  expect_warning(ggplot2::ggplot_build(overlay), NA)
  expect_warning(ggplot2::ggplot_build(checker), NA)
  expect_warning(ggplot2::ggplot_build(edge), NA)
})

test_that("standalone sagittal orientation annotations use P/A horizontally", {
  layers <- annotate_orientation("sagittal", dims = c(8L, 12L))
  labels <- vapply(
    layers,
    function(layer) as.character(layer$geom_params$grob$label),
    character(1)
  )
  expect_equal(labels, c("P", "A", "S", "I"))
})
