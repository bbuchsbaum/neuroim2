# Shared test helpers for neuroim2
# Automatically loaded by testthat before tests run.

#' Create a NeuroSpace with sensible defaults
#' @note NeuroSpace uses a 3D affine internally. For >3D spaces, only the
#'   first 3 elements of spacing/origin define the spatial transform.
make_space <- function(dim = c(10L, 10L, 10L), spacing = c(1, 1, 1),
                       origin = c(0, 0, 0), trans = NULL) {
  dim <- as.integer(dim)
  if (!is.null(trans)) {
    NeuroSpace(dim, trans = trans)
  } else {
    NeuroSpace(dim, spacing = spacing, origin = origin)
  }
}

#' Create a DenseNeuroVol with random data
make_vol <- function(dim = c(10L, 10L, 10L), data = NULL, label = "", ...) {
  sp <- make_space(dim = dim, ...)
  if (is.null(data)) data <- rnorm(prod(dim))
  DenseNeuroVol(data, sp, label = label)
}

#' Create a DenseNeuroVec (4D) with random data
make_vec <- function(dim = c(10L, 10L, 10L), ntime = 20L, data = NULL, ...) {
  sdim <- as.integer(dim)
  ntime <- as.integer(ntime)
  sp4 <- make_space(dim = c(sdim, ntime), ...)
  if (is.null(data)) data <- array(rnorm(prod(sdim) * ntime), c(sdim, ntime))
  DenseNeuroVec(data, sp4)
}

#' Create a LogicalNeuroVol mask
make_mask <- function(dim = c(10L, 10L, 10L), frac = 0.3, ...) {
  sp <- make_space(dim = dim, ...)
  m <- array(runif(prod(dim)) < frac, as.integer(dim))
  LogicalNeuroVol(m, sp)
}

#' Create a SparseNeuroVol with random nonzero entries
make_sparse_vol <- function(dim = c(10L, 10L, 10L), frac = 0.3, ...) {
  sp <- make_space(dim = dim, ...)
  n <- prod(dim)
  idx <- sort(sample.int(n, floor(n * frac)))
  sv <- Matrix::sparseVector(x = rnorm(length(idx)), i = idx, length = n)
  SparseNeuroVol(sv, sp, indices = idx)
}

#' Create a SparseNeuroVec with random data
make_sparse_vec <- function(dim = c(10L, 10L, 10L), ntime = 20L, frac = 0.3, ...) {
  sdim <- as.integer(dim)
  ntime <- as.integer(ntime)
  sp4 <- make_space(dim = c(sdim, ntime), ...)
  mask <- array(runif(prod(sdim)) < frac, sdim)
  nk <- sum(mask)
  dat <- matrix(rnorm(ntime * nk), nrow = ntime, ncol = nk)
  SparseNeuroVec(dat, sp4, mask = mask)
}

#' Write a temporary NIfTI file from a volume, return the path
make_temp_nifti <- function(vol, ext = ".nii.gz") {
  tmp <- tempfile(fileext = ext)
  write_vol(vol, tmp)
  tmp
}

# --- test helpers for version- and environment-dependent plotting -------------

#' Read a ggplot's title across ggplot2 versions.
#'
#' ggplot2 gained the exported `get_labs()` accessor in 3.5.2; before that the
#' labels are read off the plot object directly. DESCRIPTION sets no ggplot2
#' version floor, so tests must work either way.
plot_title <- function(p) {
  if ("get_labs" %in% getNamespaceExports("ggplot2")) {
    ggplot2::get_labs(p)$title
  } else {
    p$labels$title
  }
}

#' Skip visual snapshot comparisons outside the environment that produced them.
#'
#' vdiffr golden images are SVGs that encode the font metrics of the machine
#' that generated them, so byte-comparing them only succeeds on that machine.
#' The R-CMD-check matrix spans Windows, macOS and Linux, which cannot all
#' match one stored SVG, so on the other platforms these report failures that
#' say nothing about the package.
#'
#' Set NEUROIM2_SNAPSHOT_TESTS=true to run them -- do that on the reference
#' machine when reviewing intentional visual changes, and regenerate with
#' testthat::snapshot_accept() there.
skip_unless_snapshot_env <- function() {
  testthat::skip_if_not_installed("vdiffr")
  if (!identical(Sys.getenv("NEUROIM2_SNAPSHOT_TESTS"), "true")) {
    testthat::skip("visual snapshots: set NEUROIM2_SNAPSHOT_TESTS=true to compare")
  }
}
