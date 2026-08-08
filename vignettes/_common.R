# Shared setup for all neuroim2 vignettes.
#
# Sourced from the setup chunk of every article so that knitr options, the
# plotting theme, and the demonstration data are defined in exactly one place.
#
# Demonstration data policy
# -------------------------
# Three rules, applied everywhere:
#
#   * anatomy / intensity / filtering / plotting -> demo_anatomy()
#   * masks, coordinates and spatial support     -> demo_mask()
#   * anything with a time axis                  -> demo_bold()
#
# Nothing with a time axis is ever taken from the shipped `global_mask_v4.nii`:
# that file is a binary mask copied across four timepoints, so every "time
# series" drawn from it is a constant vector. demo_bold() uses simulate_fmri(),
# which gives voxels genuine AR(1) temporal structure, spatial smoothness and
# latent components.
#
# Caveat on the anatomical image: mni_downsampled.nii.gz has an inconsistent
# header. Its pixdim records 4.02 mm voxels but its sform -- which is what
# neuroim2 follows, correctly, because sform_code > 0 -- is the identity with a
# translation, so spacing(demo_anatomy()) is 1 mm and its world coordinates do
# not describe a real brain. The image is therefore used for display and
# intensity work only, where the affine affects nothing but aspect ratio (and
# it is isotropic either way). Every millimetre-coordinate example uses
# demo_mask(), whose affine is consistent: 3.5 x 3.5 x 3.7 mm, LAS.

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 4.2,
  dpi = 110,
  fig.alt = "Plot output"
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  if (requireNamespace("albersdown", quietly = TRUE)) {
    ggplot2::theme_set(albersdown::theme_albers(family = "red", preset = "homage"))
  } else {
    ggplot2::theme_set(neuroim2::theme_neuro())
  }
}

# A dark-to-light ramp that suits the shipped anatomical image.
anatomy_cmap <- c("#242424", "#9a9a9a", "#f4f4f4")

# --- demonstration data ------------------------------------------------------

#' Real anatomical image, 48 x 57 x 48, continuous intensities. Display only:
#' see the header caveat above.
demo_anatomy <- function() {
  neuroim2::read_vol(
    system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2")
  )
}

#' Binary brain mask on a 64 x 64 x 25 EPI grid (29,532 in-mask voxels), with a
#' consistent 3.5 x 3.5 x 3.7 mm affine.
demo_mask <- function() {
  neuroim2::read_vol(
    system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")
  )
}

#' Simulated BOLD series on the mask grid: AR(1) noise, 6 mm spatial
#' smoothness, latent components. ~3 seconds to generate.
demo_bold <- function(n_time = 60, seed = 1) {
  neuroim2::simulate_fmri(demo_mask(), n_time = n_time, seed = seed)
}

#' A boxcar task regressor for `n_time` scans: `block` on, `block` off.
demo_design <- function(n_time = 60, block = 10) {
  rep(rep(c(0, 1), each = block), length.out = n_time)
}

#' A brain-shaped mask on the anatomical image's own grid, so that a statistical
#' map and the anatomy it is drawn on share a NeuroSpace.
demo_anatomy_mask <- function(anat = demo_anatomy()) {
  arr <- as.array(anat)
  neuroim2::as.mask(neuroim2::NeuroVol(
    array(arr > stats::quantile(arr[arr > 0], 0.15), dim(anat)),
    neuroim2::space(anat)
  ))
}

#' A signed t-statistic map on the grid of `mask`.
#'
#' Simulates a series, plants task signal in two spheres, then actually fits the
#' design at every voxel. The result has the properties a real statistic map has:
#' signed, spatially smooth, noisy, with a null distribution of roughly unit
#' spread -- none of which a hand-drawn blob image would have.
demo_stat_map <- function(mask = demo_mask(), n_time = 30, seed = 2,
                          effect = c(1.3, -1.2), radius_mm = 12) {
  # Express smoothness in voxels rather than millimetres so the cost does not
  # explode on fine grids.
  fwhm <- 2.5 * min(neuroim2::spacing(mask))
  bold <- neuroim2::simulate_fmri(mask,
    n_time = n_time, seed = seed,
    spatial_fwhm = fwhm, hetero_fwhm = 4 * fwhm, factor_fwhm = 2 * fwhm
  )
  design <- demo_design(n_time)

  # Two centres placed either side of the mask's centroid, so this works on any
  # grid without hand-tuned voxel coordinates.
  inmask <- which(as.vector(mask) > 0)
  g <- neuroim2::index_to_grid(mask, inmask)
  mid <- round(colMeans(g))
  span <- round(0.22 * (apply(g, 2, max) - apply(g, 2, min)))
  centres <- list(mid + c(span[1], 0, 0), mid - c(span[1], 0, 0))

  Y <- as.matrix(bold)
  for (i in seq_along(centres)) {
    roi <- neuroim2::spherical_roi(mask, centres[[i]],
      radius = radius_mm, nonzero = TRUE
    )
    idx <- neuroim2::indices(roi)
    Y[idx, ] <- Y[idx, ] + effect[i] * rep(design, each = length(idx))
  }

  X <- cbind(1, design)
  XtXi <- solve(crossprod(X))
  B <- Y %*% X %*% XtXi
  s2 <- rowSums((Y - tcrossprod(B, X))^2) / (ncol(Y) - ncol(X))
  se <- sqrt(s2 * XtXi[2, 2])
  tstat <- ifelse(se > 0, B[, 2] / se, 0)
  tstat[as.vector(mask) == 0] <- 0

  neuroim2::NeuroVol(array(tstat, dim(mask)), neuroim2::space(mask))
}
