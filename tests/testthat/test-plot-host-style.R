# A host application restyles the plot family through options, without
# touching each call; unset options leave the styles as they are.

test_that("style tokens are unchanged without host options", {
  withr::local_options(neuro.ground = NULL, neuro.card = NULL, neuro.cross = NULL)
  tk <- neuro_style_tokens("report")
  expect_equal(tk$tile, "#000000")
  expect_equal(tk$card, "#f7f7f5")
  expect_equal(tk$cross, "#5fd0ff")
  expect_equal(tile_colour("grays"), resolve_cmap("grays", 2L)[[1L]])
})

test_that("host options set the ground, card and crosshair, with readable text", {
  withr::local_options(neuro.ground = "#ebe9e3", neuro.card = "#ffffff", neuro.cross = "#18181b")
  tk <- neuro_style_tokens("dark")
  expect_equal(tk$tile, "#ebe9e3")
  expect_equal(tk$tile_fg, "#1a1a1a")      # dark text on a light ground
  expect_equal(tk$card, "#ffffff")
  expect_equal(tk$fg, "#161616")
  expect_equal(tk$cross, "#18181b")
  expect_equal(tile_colour("grays"), "#ebe9e3")
})

test_that("on a host ground the background's darkest voxels are left empty", {
  withr::local_options(neuro.ground = "#ebe9e3")
  bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  p <- plot_montage(bg, zlevels = 24, colorbar = FALSE)
  expect_s3_class(p, "neuro_fig")
  withr::local_options(neuro.ground = NULL)
  expect_s3_class(plot_montage(bg, zlevels = 24, colorbar = FALSE), "neuro_fig")
})
