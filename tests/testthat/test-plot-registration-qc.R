library(testthat)

make_registration_qc_volumes <- function(dims = c(8L, 9L, 5L), space = neuroim2::NeuroSpace(dims)) {
  list(
    bg = neuroim2::NeuroVol(array(seq_len(prod(dims)), dim = dims), space),
    ov = neuroim2::NeuroVol(array(rev(seq_len(prod(dims))), dim = dims), space),
    e1 = neuroim2::NeuroVol(array(seq_len(prod(dims)) %% 3L, dim = dims), space),
    e2 = neuroim2::NeuroVol(array(seq_len(prod(dims)) %% 5L, dim = dims), space)
  )
}

test_that("plot_checkerboard returns ggplot panels and an assembled figure", {
  vols <- make_registration_qc_volumes()

  result <- neuroim2::plot_checkerboard(
    vols$bg,
    vols$ov,
    zlevels = c(2L, 4L),
    tile = 2L,
    ncol = 2L,
    title = "Checker",
    draw = FALSE,
    assemble = FALSE
  )

  expect_length(result, 2L)
  expect_true(all(vapply(result, inherits, logical(1), what = "ggplot")))
  # Panels are labelled with the world coordinate of the slice: voxel 2 on a
  # unit grid with zero origin sits at z = 1 mm.
  expect_world_label(tile_slice_label(result[[1L]]), "z = 1")
  expect_world_label(tile_slice_label(result[[2L]]), "z = 3")
  expect_equal(attr(result, "labels")$title, "Checker")

  # The default return is one assembled figure: visible when not drawn,
  # invisible after draw = TRUE prints it.
  fig <- neuroim2::plot_checkerboard(vols$bg, vols$ov, zlevels = c(2L, 4L),
                                     tile = 2L, draw = FALSE)
  expect_s3_class(fig, "patchwork")
  expect_visible(neuroim2::plot_checkerboard(vols$bg, vols$ov, zlevels = 2L,
                                             tile = 2L, draw = FALSE))
  grDevices::pdf(NULL)
  expect_invisible(neuroim2::plot_checkerboard(vols$bg, vols$ov, zlevels = 2L,
                                               tile = 2L, draw = TRUE))
  grDevices::dev.off()
})

test_that("plot_edge_overlay returns ggplot panels and an assembled figure", {
  vols <- make_registration_qc_volumes()

  result <- neuroim2::plot_edge_overlay(
    vols$bg,
    vols$e1,
    vols$e2,
    zlevels = c(2L, 4L),
    edge_thresh = 0,
    ncol = 2L,
    title = "Edges",
    draw = FALSE,
    assemble = FALSE
  )

  expect_length(result, 2L)
  expect_true(all(vapply(result, inherits, logical(1), what = "ggplot")))
  expect_world_label(tile_slice_label(result[[1L]]), "z = 1")
  expect_world_label(tile_slice_label(result[[2L]]), "z = 3")
  expect_equal(attr(result, "labels")$title, "Edges")

  fig <- neuroim2::plot_edge_overlay(vols$bg, vols$e1, vols$e2,
                                     zlevels = c(2L, 4L), draw = FALSE)
  expect_s3_class(fig, "patchwork")
  expect_visible(neuroim2::plot_edge_overlay(vols$bg, vols$e1, vols$e2,
                                             zlevels = 2L, draw = FALSE))
  grDevices::pdf(NULL)
  expect_invisible(neuroim2::plot_edge_overlay(vols$bg, vols$e1, vols$e2,
                                               zlevels = 2L, draw = TRUE))
  grDevices::dev.off()
})

test_that("registration QC plots reject volumes on different grids", {
  dims <- c(8L, 9L, 5L)
  vols <- make_registration_qc_volumes(dims)
  shifted_space <- neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2))
  shifted <- neuroim2::NeuroVol(array(0, dim = dims), shifted_space)

  expect_error(
    neuroim2::plot_checkerboard(vols$bg, shifted, zlevels = 2L, draw = FALSE),
    "same NeuroSpace grid"
  )
  expect_error(
    neuroim2::plot_edge_overlay(vols$bg, vols$e1, shifted, zlevels = 2L, draw = FALSE),
    "same NeuroSpace grid"
  )
})

test_that("registration QC plots validate panel layout arguments", {
  vols <- make_registration_qc_volumes()

  expect_error(
    neuroim2::plot_checkerboard(vols$bg, vols$ov, zlevels = integer(0), draw = FALSE),
    "`zlevels`"
  )
  expect_error(
    neuroim2::plot_checkerboard(vols$bg, vols$ov, zlevels = 2L, ncol = 0L),
    "`ncol`"
  )
  expect_error(
    neuroim2::plot_edge_overlay(vols$bg, vols$e1, vols$e2, zlevels = 99L, draw = FALSE),
    "`zlevels`"
  )
  expect_error(
    neuroim2::plot_edge_overlay(vols$bg, vols$e1, vols$e2, zlevels = 2L, ncol = 0L),
    "`ncol`"
  )
})

test_that("registration QC plots draw panel grids with layout labels", {
  vols <- make_registration_qc_volumes(c(5L, 6L, 4L))
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit(unlink(tf), add = TRUE)

  expect_error(
    neuroim2::plot_checkerboard(
      vols$bg,
      vols$ov,
      zlevels = c(1L, 3L),
      ncol = 2L,
      title = "Checker",
      subtitle = "Sub",
      caption = "Cap",
      draw = TRUE
    ),
    NA
  )
  expect_error(
    neuroim2::plot_edge_overlay(
      vols$bg,
      vols$e1,
      vols$e2,
      zlevels = c(1L, 3L),
      ncol = 2L,
      title = "Edges",
      subtitle = "Sub",
      caption = "Cap",
      draw = TRUE
    ),
    NA
  )
  grDevices::dev.off()
  expect_gt(file.info(tf)$size, 0)
})

test_that("plot_edge_overlay keeps all-zero edge slices transparent", {
  dims <- c(4L, 4L, 2L)
  sp <- neuroim2::NeuroSpace(dims)
  bg <- neuroim2::NeuroVol(array(1, dim = dims), sp)
  edges <- neuroim2::NeuroVol(array(0, dim = dims), sp)

  p <- neuroim2::plot_edge_overlay(bg, edges, edges, zlevels = 1L, draw = FALSE,
                                   assemble = FALSE)[[1L]]
  # The fixed, moving and agreement channels are annotation_custom raster
  # layers over the background; with no edges all must be transparent.
  edge_layers <- Filter(function(l) layer_geom(l) == "GeomCustomAnn", p$layers)
  expect_length(edge_layers, 3L)
  for (layer in edge_layers) {
    edge_grob <- layer$geom_params$grob
    alpha_hex <- unique(substr(as.character(edge_grob$raster), 8L, 9L))
    expect_equal(alpha_hex, "00")
  }
})
