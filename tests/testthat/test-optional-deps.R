library(testthat)
library(neuroim2)

# mmap and bigstatsr are Suggests: neuroim2 must work without them (e.g. webR).

test_that("masked reads fall back to streaming when mmap is unavailable", {
  fname <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  mask <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")) > 0
  ref <- read_vec(fname, mask = mask)

  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  streamed <- read_vec(fname, mask = mask)
  expect_s4_class(streamed, "SparseNeuroVec")
  expect_equal(as.matrix(streamed), as.matrix(ref))
})

test_that("memory-mapped modes fail informatively when mmap is unavailable", {
  fname <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  expect_error(read_vec(fname, mode = "mmap"), "Package 'mmap' is required")
  expect_error(MappedNeuroVec(fname), "Package 'mmap' is required")
})

test_that("BigNeuroVec validity requires an FBM in the data slot", {
  sp <- NeuroSpace(c(2, 2, 1, 2))
  mask <- LogicalNeuroVol(array(TRUE, c(2, 2, 1)), drop_dim(sp))
  expect_error(
    new("BigNeuroVec", space = sp, mask = mask,
        map = IndexLookupVol(space(mask), 1:4), data = matrix(0, 2, 4)),
    "FBM"
  )
})

test_that("as_mmap stops before writing a file when mmap is unavailable", {
  vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
  out <- tempfile(fileext = ".nii")
  local_mocked_bindings(.has_mmap = function() FALSE, .package = "neuroim2")
  expect_error(as_mmap(vec, file = out), "Package 'mmap' is required")
  expect_false(file.exists(out))
  sparse <- as.sparse(vec, vec[[1]] > 0)
  expect_error(as_mmap(sparse, file = out), "Package 'mmap' is required")
  expect_false(file.exists(out))
})

test_that("assemble = FALSE returns plain ggplot panels without patchwork", {
  bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  stat <- bg
  stat[] <- (bg[] > stats::quantile(bg[], 0.7)) * 3
  moving <- bg
  local_mocked_bindings(
    .require_patchwork = function() stop("patchwork unavailable (mocked)"),
    .package = "neuroim2"
  )
  is_panels <- function(p) is.list(p) && length(p) > 0 &&
    all(vapply(p, inherits, logical(1), "ggplot"))
  expect_true(is_panels(plot_ortho(bg, assemble = FALSE)))
  expect_true(is_panels(plot_overlay(bg, stat, ov_thresh = 1, assemble = FALSE)))
  expect_true(is_panels(plot_checkerboard(bg, moving, assemble = FALSE)))
  expect_true(is_panels(plot_edge_overlay(bg, bg, moving, compute_edges = TRUE,
                                          assemble = FALSE)))
  expect_error(plot_ortho(bg), "patchwork unavailable")
})
