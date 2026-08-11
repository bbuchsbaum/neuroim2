#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Spatial Filtering Methods for Neuroimaging Data
#'
#' @name spatial-filter
#' @description Methods for applying spatial filters to neuroimaging data
NULL

#' @importFrom methods new
#' @importFrom stats dnorm
NULL

# ---------------------------------------------------------------------------
# Gaussian blur engine choice
#
# Two implementations, identical up to floating-point accumulation order:
#   * gaussian_blur_cpp()     -- the full (2w+1)^3 kernel, evaluated only at mask
#                                voxels. No working memory.
#   * gaussian_blur_sep_cpp() -- three 1-D passes over the whole volume, which is
#                                3(2w+1) taps instead of (2w+1)^3, but touches
#                                every voxel and needs 2-3 scratch vectors of
#                                prod(dim) doubles.
#
# Which is cheaper depends on the window *and* the mask fraction, so estimate
# both rather than gating on the window alone. Tap counts alone give
#
#   dense       ~ (2w+1)^3 * n_mask
#   separable   ~ passes * (2w+1) * n_vox     (passes = 3, or 6 when normalising,
#                                              which needs blur(y) and blur(m))
#
# which would put the crossover at  n_mask/n_vox  =  passes / (2w+1)^2. That is
# far too conservative in practice, because the two kernels do not cost the same
# per tap: the separable passes walk contiguous memory with no per-tap branch,
# while the dense kernel pays a bounds test and a mask lookup on every tap, and
# that gap widens with the window. Measuring the actual crossover on a 91x109x91
# volume over mask fractions 0.02-1.00 (both contiguous and scattered masks) put
# it at roughly
#
#   normalize = TRUE   w = 1: 0.15   w = 2: 0.09   w = 3: 0.05   w = 4: 0.04
#   normalize = FALSE  w = 1: 0.08   w = 2: 0.05   w = 3: 0.04   w = 4: 0.02
#
# i.e. ~ 0.15 * passes / (2w+1)^1.5. Over that 224-cell grid this rule picks the
# slower kernel by more than 5% in 7 cells and never by more than 1.40x (on a
# call costing 14ms), while its total time is within 0.6% of always choosing the
# per-cell winner -- against 28% for the tap-count rule, which gives away
# speedups of up to 5x on ordinary brain-mask-sized inputs, and 8.8x for always
# taking the dense path.
#
# The constant is a calibration, not a law: it is a ratio of per-tap costs, so it
# drifts with the compiler and the machine. An independent measurement on other
# hardware put the w = 1, normalize = TRUE crossover nearer 0.20 than the 0.15
# above, against a threshold of 0.173 -- so a typical brain mask can land in a
# band where this fires ~1.2x early. That spread is the honest width of the fit,
# and it is tolerable because near the crossover the two kernels cost nearly the
# same by construction. The thresholds are pinned in
# tests/testthat/test-gaussian-blur-separable.R so re-calibrating is deliberate.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.gaussian_blur_prefers_separable <- function(window, n_mask, n_vox, normalize) {
  # Total in its arguments: this only chooses an engine, so a degenerate window
  # or an empty volume must fall through to the kernel's own validation rather
  # than error here.
  if (!isTRUE(n_vox > 0) || length(window) != 1L ||
      !isTRUE(is.finite(window)) || !isTRUE(window >= 1)) {
    return(FALSE)
  }
  passes <- if (isTRUE(normalize)) 6 else 3
  (n_mask / n_vox) > (0.15 * passes / (2 * window + 1)^1.5)
}

#' Gaussian blur, dispatching to whichever kernel is cheaper
#'
#' @keywords internal
#' @noRd
.gaussian_blur_engine <- function(arr, mask_idx, window, sigma, spacing,
                                  normalize = TRUE, full_mask = FALSE) {
  window <- as.integer(window)
  n_vox <- length(arr)
  # `full_mask = TRUE` means "every voxel", and lets the caller skip building an
  # n_vox-long index vector it does not otherwise need.
  n_mask <- if (isTRUE(full_mask)) n_vox else length(mask_idx)
  mask_idx <- if (isTRUE(full_mask)) integer(0) else as.integer(mask_idx)

  # The separable kernel requires exactly three dimensions; the dense one blurs
  # the first volume of a higher-dimensional array. Nothing reachable from an
  # exported function passes one, but the engine must not turn that into an
  # error where it used to be a result.
  is_3d <- length(dim(arr)) == 3L

  if (is_3d &&
      .gaussian_blur_prefers_separable(window, n_mask, n_vox, normalize)) {
    gaussian_blur_sep_cpp(arr, mask_idx, window, sigma, spacing, normalize,
                          isTRUE(full_mask))
  } else {
    # The dense kernel has no full-mask convention, so hand it the indices.
    if (isTRUE(full_mask)) mask_idx <- seq_len(n_vox)
    gaussian_blur_cpp(arr, mask_idx, window, sigma, spacing, normalize)
  }
}

#' Kernel half-width, in voxels, that keeps `truncate` standard deviations
#'
#' `sigma` is in millimetres and `window` counts voxels, so the conversion needs
#' the voxel size. The largest value over the three axes is used rather than a
#' per-axis width: the kernel builders take one window, and the extra taps a
#' coarse axis then evaluates lie beyond `truncate` sigma and carry negligible
#' weight. The rounding is `scipy.ndimage.gaussian_filter`'s -- it takes
#' `int(truncate * sd + 0.5)` voxels, with its default `truncate = 4.0` -- so
#' the kernel matches the reference implementation tap for tap on isotropic
#' data rather than overshooting it.
#' @keywords internal
#' @noRd
.gaussian_window <- function(sigma, spacing, truncate = 4) {
  sp <- as.numeric(spacing)
  sp <- sp[is.finite(sp) & sp > 0]
  if (!length(sp)) sp <- 1
  max(1L, as.integer(floor(truncate * sigma / min(sp) + 0.5)))
}

#' Gaussian Blur for Volumetric Images
#'
#' @description
#' This function applies an isotropic discrete Gaussian kernel to smooth a volumetric image (3D brain MRI data).
#' The blurring is performed within a specified image mask, with customizable kernel parameters.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask.
#'   This mask defines the region where the blurring is applied. If not provided, the entire volume is processed.
#'   It must have the same dimensions as \code{vol}.
#' @param sigma A single finite positive number specifying the standard deviation of the Gaussian
#'   kernel in the \strong{spatial units of the image} (millimetres for a typical NIfTI), not
#'   voxels. The kernel is evaluated in physical distance using \code{spacing(vol)}, so the same
#'   \code{sigma} produces the same physical smoothing regardless of voxel size. Default is 2;
#'   ignored when \code{fwhm} is given.
#' @param fwhm Full width at half maximum of the kernel in the image's spatial units, the form in
#'   which smoothing is usually specified. Supplying it sets \code{sigma} to
#'   \code{fwhm / (2 * sqrt(2 * log(2)))}; giving both is an error.
#' @param window Kernel half-width in \emph{voxels}: the number of voxels included on each side
#'   of the centre, so \code{window = 1} is a 3x3x3 kernel. \code{NULL} (the default) derives it
#'   from \code{sigma}, \code{truncate}, and the voxel spacing, using the largest half-width over
#'   the three axes so no axis is clipped early. Supplying a number instead truncates the kernel
#'   there, which can under-smooth --- see \emph{Details}. A window larger than the volume is
#'   clamped to \code{max(dim(vol))}, which changes no result: every kernel tap that far out is
#'   out of bounds from every voxel.
#' @param truncate How many standard deviations of the Gaussian to keep when \code{window} is
#'   derived. Default 4, matching \code{scipy.ndimage.gaussian_filter}. Smaller values are faster
#'   and smooth less than requested.
#' @param normalize A logical value controlling how the mask boundary is handled. When \code{TRUE}
#'   (the default), the blur is \emph{insulated} to the mask: each in-mask output voxel is computed
#'   from in-mask neighbours only and the kernel is renormalized by the in-mask weight (a
#'   "smooth-in-mask" convolution, cf. AFNI \code{3dBlurInMask}). Out-of-mask values --- finite or
#'   not (including \code{NaN}/\code{Inf}) --- never influence in-mask outputs, and edge voxels are
#'   not biased toward the mask exterior. When \code{FALSE}, the legacy behavior is used: the full
#'   kernel is applied with zero padding outside the volume and out-of-mask neighbour values are
#'   read into the convolution (retained only for backward compatibility).
#'
#' @return A \code{\linkS4class{NeuroVol}} object representing the smoothed image.
#'
#' @details
#' The function uses a C++ implementation for efficient Gaussian blurring. The blurring is applied
#' only to voxels within the specified mask (or the entire volume if no mask is provided).
#'
#' \strong{Kernel width.} A Gaussian has infinite support, so it has to be cut off somewhere, and
#' where it is cut off changes the smoothing that is actually applied. \code{window} is measured
#' in voxels while \code{sigma} is in millimetres, so a fixed \code{window} truncates at a
#' different point on every acquisition. The old default of \code{window = 1} cut a
#' \code{sigma = 2} mm kernel off at one standard deviation on 2 mm voxels and delivered 3.5 mm
#' FWHM where 4.7 mm was asked for --- and 1.9 mm on 1 mm voxels. Leaving \code{window = NULL}
#' derives it from \code{sigma} instead, so the smoothing you ask for is the smoothing you get:
#'
#' \tabular{lrrr}{
#'   \strong{call} \tab \strong{delivered} \tab \strong{requested} \tab \strong{shortfall} \cr
#'   \code{sigma = 2, window = 1} \tab 3.49 mm \tab 4.71 mm \tab 26\% \cr
#'   \code{sigma = 4, window = 1} \tab 3.76 mm \tab 9.42 mm \tab 60\% \cr
#'   \code{sigma = 4} (derived)   \tab 9.42 mm \tab 9.42 mm \tab 0\% \cr
#' }
#'
#' Passing \code{window} explicitly still truncates exactly as before, so old calls are
#' reproducible; it is only the default that changed.
#'
#' \strong{Units.} \code{sigma} is in the spatial units of the image (millimetres for typical
#' NIfTI data), \emph{not} voxels: neighbour offsets are multiplied by \code{spacing(vol)} before
#' the kernel is evaluated. \code{window}, by contrast, is a half-width counted in voxels. Mixing
#' the two is the easy mistake --- converting an FWHM to sigma in voxels and passing that as
#' \code{sigma} silently under-smooths, with no error or warning. To smooth to a target FWHM:
#'
#' \preformatted{
#' fwhm_mm  <- 6
#' sigma_mm <- fwhm_mm / (2 * sqrt(2 * log(2)))          # ~2.55 mm
#' window   <- ceiling(3 * sigma_mm / mean(spacing(vol))) # cover +/- 3 sigma, in voxels
#' gaussian_blur(vol, mask, sigma = sigma_mm, window = window)
#' }
#'
#' \code{window} must be large enough to contain the kernel: because the kernel is truncated at
#' \code{window} voxels, a \code{window} much smaller than \code{3 * sigma_mm / spacing} clips the
#' Gaussian tails and delivers less smoothing than \code{sigma} implies.
#'
#' With the default \code{normalize = TRUE}, smoothing a masked statistical map (e.g. a first-level
#' fMRI coefficient map written as \code{NaN} outside the brain) within its brain mask preserves the
#' full in-mask coverage; out-of-mask \code{NaN}s do not erode the masked region and the renormalized
#' edge values match the interior on their shared support.
#'
#' @examples
#' # Load a sample brain mask
#' brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Apply Gaussian blurring to the brain volume
#' blurred_vol <- gaussian_blur(brain_mask, brain_mask, sigma = 2, window = 1)
#'
#' # View a slice of the original and blurred volumes
#' image(brain_mask[,,12])
#' image(blurred_vol[,,12])
#'
#' # Smooth to a target FWHM. `sigma` is in millimetres, so convert FWHM -> sigma
#' # in mm (NOT voxels); `window` is a voxel half-width sized to cover +/- 3 sigma.
#' fwhm_mm  <- 6
#' sigma_mm <- fwhm_mm / (2 * sqrt(2 * log(2)))
#' window   <- ceiling(3 * sigma_mm / mean(spacing(brain_mask)))
#' smoothed <- gaussian_blur(brain_mask, brain_mask, sigma = sigma_mm, window = window)
#'
#' @seealso
#' \code{\link{NeuroVol-class}}, \code{\link{LogicalNeuroVol-class}}, \code{\link{bilateral_filter}}
#'
#' @references
#' Gaussian blur: https://en.wikipedia.org/wiki/Gaussian_blur
#'
#' @export
gaussian_blur <- function(vol, mask, sigma = 2, window = NULL, normalize = TRUE,
                          fwhm = NULL, truncate = 4) {
  if (!inherits(vol, "NeuroVol") && !inherits(vol, "NeuroVec")) {
    cli::cli_abort("{.arg vol} must be a {.cls NeuroVol} or {.cls NeuroVec} object.")
  }
  if (!is.null(fwhm)) {
    if (!missing(sigma)) {
      cli::cli_abort(c("Give either {.arg sigma} or {.arg fwhm}, not both.",
                       "i" = "{.arg fwhm} is {.code sigma * 2 * sqrt(2 * log(2))}."))
    }
    if (!is.numeric(fwhm) || length(fwhm) != 1L || !is.finite(fwhm) || fwhm <= 0) {
      cli::cli_abort("{.arg fwhm} must be a single finite positive number, not {.val {fwhm}}.")
    }
    sigma <- fwhm / (2 * sqrt(2 * log(2)))
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    cli::cli_abort("{.arg sigma} must be a single finite positive number, not {.val {sigma}}.")
  }
  if (!is.numeric(truncate) || length(truncate) != 1L || !is.finite(truncate) ||
      truncate <= 0) {
    cli::cli_abort("{.arg truncate} must be a single finite positive number, not {.val {truncate}}.")
  }
  if (is.null(window)) {
    window <- .gaussian_window(sigma, spacing(vol)[1:3], truncate)
  }
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window < 1) {
    cli::cli_abort("{.arg window} must be a single finite number >= 1, not {.val {window}}.")
  }
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    cli::cli_abort("{.arg normalize} must be a single {.cls logical} (TRUE or FALSE).")
  }
  spatial_dim <- as.integer(dim(vol)[1:3])
  if (!missing(mask)) {
    if (!inherits(mask, "NeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls NeuroVol} object.")
    }
    # A 4-D image is smoothed volume by volume, so the mask is 3-D either way.
    if (!identical(as.integer(dim(mask)), spatial_dim)) {
      cli::cli_abort(c(
        "{.arg mask} and {.arg vol} must have the same spatial dimensions.",
        i = "{.arg vol} is {.val {dim(vol)}}, {.arg mask} is {.val {dim(mask)}}."
      ))
    }
  }

  # Every kernel tap at least max(dim) voxels away is out of bounds from every
  # voxel, so clipping the window here changes no result. It does keep absurd
  # windows out of the kernel builders, where (2w+1)^3 used to overflow.
  # min() first: as.integer(1e10) is NA, not a clamped value.
  window <- as.integer(min(window, max(spatial_dim)))

  nvox <- prod(spatial_dim)
  full_mask <- missing(mask)
  mask.idx <- if (full_mask) integer(0) else which(mask != 0)
  vox_spacing <- spacing(vol)[1:3]

  if (inherits(vol, "NeuroVol")) {
    target_space <- if (missing(mask)) space(vol) else space(mask)
    farr <- .gaussian_blur_engine(.image_array(vol), mask.idx, window, sigma,
                                  vox_spacing, normalize, full_mask)
    return(NeuroVol(farr, target_space))
  }

  # 4-D: smooth each volume with the same kernel. Which engine is cheaper does
  # not depend on the volume, so it is decided once here rather than re-derived
  # on every iteration.
  nvols <- dim(vol)[4]
  n_mask <- if (full_mask) nvox else length(mask.idx)
  if (.gaussian_blur_prefers_separable(window, n_mask, nvox, normalize)) {
    # The compiled driver walks the run itself: no volume is copied out and no
    # result copied back, the scratch buffers are allocated per worker instead
    # of per volume, and blur(m) -- which is the same for every volume -- is
    # computed once. It also runs volumes in parallel.
    arr <- .image_array(vol)
    if (length(dim(arr)) != 4L) dim(arr) <- c(spatial_dim, nvols)
    out <- gaussian_blur_sep_4d_cpp(arr, if (full_mask) integer(0) else as.integer(mask.idx),
                                    window, sigma, vox_spacing, normalize, full_mask)
    return(DenseNeuroVec(out, space(vol), label = vol@label,
                         volume_labels = volume_labels(vol)))
  }

  # The dense kernel is chosen when the mask is a small fraction of the volume.
  # It has no 4-D form, so the loop stays -- but everything invariant across
  # volumes is still hoisted out of it. Doing this in user code instead --
  # gaussian_blur(x[[i]], ...) in a loop -- would re-run the argument checks,
  # rebuild a NeuroVol and re-derive the mask once per volume.
  flat <- .image_array(vol)
  dim(flat) <- c(nvox, nvols)
  out <- numeric(nvox * nvols)
  dim(out) <- c(nvox, nvols)
  buf_dim <- spatial_dim
  for (i in seq_len(nvols)) {
    buf <- flat[, i]
    dim(buf) <- buf_dim
    out[, i] <- .gaussian_blur_engine(buf, mask.idx, window, sigma,
                                      vox_spacing, normalize, full_mask)
  }
  dim(out) <- c(spatial_dim, nvols)
  DenseNeuroVec(out, space(vol), label = vol@label,
                volume_labels = volume_labels(vol))
}

#' Edge-Preserving Guided Filter for Volumetric Images
#'
#' @description
#' This function applies a guided filter to a volumetric image (3D brain MRI data)
#' to perform edge-preserving smoothing. The guided filter smooths the image while
#' preserving edges, providing a balance between noise reduction and structural preservation.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be filtered.
#' @param radius An integer specifying the spatial radius of the filter. Default is 4.
#' @param epsilon A numeric value specifying the regularization parameter. It controls
#'   the degree of smoothing and edge preservation. Default is 0.49 (0.7^2).
#'
#' @return A \code{\linkS4class{NeuroVol}} object representing the filtered image.
#'
#' @details
#' The guided filter operates by computing local linear models between the guidance
#' image (which is the same as the input image in this implementation) and the output.
#' The 'radius' parameter determines the size of the local neighborhood, while 'epsilon'
#' controls the smoothness of the filter.
#'
#' The implementation uses box blur operations for efficiency, which approximates
#' the behavior of the original guided filter algorithm.
#'
#' @examples
#' # Load an example brain volume
#' brain_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#'
#' # Apply guided filtering to the brain volume
#' \donttest{
#' filtered_vol <- guided_filter(brain_vol, radius = 4, epsilon = 0.49)
#'
#' # Visualize a slice of the original and filtered volumes
#' oldpar <- par(mfrow = c(1, 2))
#' image(brain_vol[,,12], main = "Original")
#' image(filtered_vol[,,12], main = "Filtered")
#' par(oldpar)
#' }
#'
#' @references
#' He, K., Sun, J., & Tang, X. (2013). Guided Image Filtering. IEEE Transactions
#' on Pattern Analysis and Machine Intelligence, 35(6), 1397-1409.
#'
#' @seealso
#' \code{\link{gaussian_blur}}, \code{\link{bilateral_filter}}, \code{\link{NeuroVol-class}}
#'
#' @export
guided_filter <- function(vol, radius = 4, epsilon = 0.7^2) {
  if (!inherits(vol, "NeuroVol")) {
    cli::cli_abort("{.arg vol} must be a {.cls NeuroVol} object.")
  }
  if (radius < 1) {
    cli::cli_abort("{.arg radius} must be >= 1, not {.val {radius}}.")
  }
  if (epsilon <= 0) {
    cli::cli_abort("{.arg epsilon} must be positive, not {.val {epsilon}}.")
  }

  mask_idx <- which(vol !=0)
  mean_I = box_blur(vol, mask_idx, radius)
  mean_II = box_blur(vol*vol, mask_idx, radius)
  var_I = mean_II - mean_I * mean_I
  mean_p = box_blur(vol, mask_idx, radius)
  mean_Ip = box_blur(vol*vol, mask_idx, radius)

  cov_Ip = mean_Ip - mean_I * mean_p
  a = cov_Ip / (var_I + epsilon)
  b = mean_p - a * mean_I
  mean_a = box_blur(a, mask_idx, radius)
  mean_b = box_blur(b, mask_idx, radius)
  out = mean_a * vol + mean_b
  ovol = NeuroVol(out, space(vol))
  ovol
}

#' Apply a bilateral filter to a volumetric image
#'
#' This function smooths a volumetric image (3D brain MRI data) using a bilateral filter.
#' The bilateral filter considers both spatial closeness and intensity similarity for smoothing.
#' Only in-mask, in-bounds neighbors contribute to each local weighted average.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask that defines the region where the filtering is applied. If not provided, the entire volume is considered.
#' @param spatial_sigma A numeric value specifying the standard deviation of the spatial Gaussian kernel (default is 2).
#' @param intensity_sigma A numeric value specifying the standard deviation of the intensity Gaussian kernel (default is 1).
#' @param window An integer specifying the number of voxels around the center voxel to include on each side. For example, window=1 for a 3x3x3 kernel (default is 1).
#' @param range_scale Optional positive numeric range scale used by the intensity kernel.
#'   If \code{NULL}, the scale is estimated as the standard deviation of the current
#'   input values inside \code{mask}. Supply a
#'   fixed value to apply the same range bandwidth across observed and null maps.
#'
#' @return A smoothed image of class \code{\linkS4class{NeuroVol}}.
#'
#' @examples
#' brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Apply bilateral filtering to the brain volume
#' filtered_vol <- bilateral_filter(brain_mask, brain_mask, spatial_sigma = 2,
#' intensity_sigma = 25, window = 1)
#'
#' @export
bilateral_filter <- function(vol, mask, spatial_sigma=2, intensity_sigma=1, window=1,
                             range_scale = NULL) {
  if (window < 1) {
    cli::cli_abort("{.arg window} must be >= 1, not {.val {window}}.")
  }
  if (spatial_sigma <= 0) {
    cli::cli_abort("{.arg spatial_sigma} must be positive, not {.val {spatial_sigma}}.")
  }
  if (intensity_sigma <= 0) {
    cli::cli_abort("{.arg intensity_sigma} must be positive, not {.val {intensity_sigma}}.")
  }
  if (!is.null(range_scale)) {
    if (!is.numeric(range_scale) || length(range_scale) != 1 || !is.finite(range_scale) || range_scale <= 0) {
      cli::cli_abort("{.arg range_scale} must be {.code NULL} or a single positive finite number.")
    }
  }
  if (!missing(mask)) {
    if (!inherits(mask, "NeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls NeuroVol} object.")
    }
  }

  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vol)))
    target_space <- space(vol)
  } else {
    mask.idx <- which(mask!=0)
    target_space <- space(mask)
  }

  arr <- as.array(vol)
  farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window),
                               spatial_sigma, intensity_sigma, spacing(vol),
                               if (is.null(range_scale)) NA_real_ else range_scale)

  out <- NeuroVol(farr, target_space)
  out
}

#' Apply a bilateral filter to each volume of a NeuroVec
#'
#' This function applies a bilateral filter to each volume of a NeuroVec object.
#' The filter is applied using a specified spatial and intensity sigma, and a given window size.
#'
#' @param vec A NeuroVec object containing the volumes to be filtered.
#' @param mask A binary mask specifying the region of interest. If not provided, the whole volume is considered.
#' @param spatial_sigma The spatial sigma for the bilateral filter (default = 2).
#' @param intensity_sigma The intensity sigma for the bilateral filter (default = 1).
#' @param window The size of the window for the bilateral filter (default = 1).
#' @param range_scale Optional positive numeric range scale. If \code{NULL}, each
#'   volume estimates its own masked standard deviation.
#' @return A NeuroVec object with the filtered volumes.
#' @examples
#' brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' out <- bilateral_filter_vec(vec,brain_mask)
#' @noRd
bilateral_filter_vec <- function(vec, mask, spatial_sigma=2, intensity_sigma=1, window=1,
                                 range_scale = NULL) {
  if (!inherits(vec, "NeuroVec")) {
    cli::cli_abort("{.arg vec} must be a {.cls NeuroVec} object.")
  }
  if (window < 1) {
    cli::cli_abort("{.arg window} must be >= 1, not {.val {window}}.")
  }
  if (spatial_sigma <= 0) {
    cli::cli_abort("{.arg spatial_sigma} must be positive, not {.val {spatial_sigma}}.")
  }
  if (intensity_sigma <= 0) {
    cli::cli_abort("{.arg intensity_sigma} must be positive, not {.val {intensity_sigma}}.")
  }
  if (!is.null(range_scale)) {
    if (!is.numeric(range_scale) || length(range_scale) != 1 || !is.finite(range_scale) || range_scale <= 0) {
      cli::cli_abort("{.arg range_scale} must be {.code NULL} or a single positive finite number.")
    }
  }

  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vec)[1:3]))
    target_space <- space(vec[[1]])
  } else {
    if (!inherits(mask, "NeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls NeuroVol} object.")
    }
    mask.idx <- which(mask!=0)
    target_space <- space(mask)
  }

  res<- lapply(seq_len(dim(vec)[4]), function(i) {
    vol_i <- vec[[i]]
    arr <- as.array(vol_i)
    farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window),
                                 spatial_sigma, intensity_sigma, spacing(vec)[1:3],
                                 if (is.null(range_scale)) NA_real_ else range_scale)
    NeuroVol(farr, target_space)
  })

  do.call(concat,res)

}

#' Apply a 4D bilateral filter to a NeuroVec
#'
#' This function applies a full 4D bilateral filter to a \code{NeuroVec},
#' smoothing jointly across space (x, y, z) and time (t). The filter uses
#' spatial, temporal, and intensity kernels to preserve edges while reducing
#' noise, leveraging a parallel C++ backend for performance.
#' Only in-mask, in-bounds spatial neighbors contribute to each local weighted average.
#'
#' @param vec A \code{\linkS4class{NeuroVec}} object (4D image).
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} or \code{\linkS4class{NeuroVol}}
#'   specifying the spatial region to process. If omitted, the entire spatial
#'   extent is processed.
#' @param spatial_sigma Numeric; standard deviation of the spatial Gaussian (default 2).
#' @param intensity_sigma Numeric; standard deviation of the intensity Gaussian (default 1).
#' @param temporal_sigma Numeric; standard deviation of the temporal Gaussian (default 1).
#' @param spatial_window Integer; half-width of the spatial window in voxels (default 1),
#'   e.g., 1 => 3x3x3 spatial neighborhood.
#' @param temporal_window Integer; half-width of the temporal window in frames (default 1),
#'   e.g., 1 => 3 timepoints (t-1, t, t+1).
#' @param temporal_spacing Numeric; spacing of the temporal dimension (e.g., TR in seconds).
#'   Default is 1. This sets the temporal scale used for the temporal kernel.
#' @param range_scale Optional positive numeric range scale used by the intensity kernel.
#'   If \code{NULL}, the scale is estimated from all finite input values inside
#'   \code{mask} across time. Supply a fixed value to use the same range bandwidth
#'   across observed and null data.
#'
#' @details
#' Parameter guidance and units:
#' - spatial_sigma: Measured in physical units (millimeters). Distances are
#'   computed using \code{spacing(vec)[1:3]}, so choose \code{spatial_sigma}
#'   relative to voxel size. As a rule of thumb, set it to about 1-2 voxel sizes
#'   (e.g., 2-4 mm for 2 mm isotropic data) for moderate smoothing.
#' - intensity_sigma: Dimensionless multiplier of the global intensity standard
#'   deviation. Internally, the filter uses exp(-(dI)^2 / (2 * (intensity_sigma * sigma_I)^2)),
#'   where sigma_I is the standard deviation of all finite voxel intensities within
#'   the mask across time. Start with 1.0 for moderate smoothing; use 0.5-0.8 to
#'   preserve more edges, or 1.5-2.0 for stronger smoothing.
#' - temporal_sigma: Measured in \code{temporal_spacing} units (e.g., seconds).
#'   Typical values are 0.5-2 x TR. Larger values blend more across time.
#'
#' Choosing the neighborhood window sizes:
#' - spatial_window controls the discrete spatial support. A common choice is
#'   \code{ceiling(2 * spatial_sigma / min(spacing(vec)[1:3]))}, which covers
#'   ~95% of a Gaussian's mass.
#' - temporal_window similarly can be set to \code{ceiling(2 * temporal_sigma / temporal_spacing)}.
#'
#' Quick presets (typical fMRI with 2-3 mm voxels and TR~2s):
#' - Light: spatial_sigma = 1 x min(spacing), intensity_sigma = 0.8,
#'   temporal_sigma = 0.5 x TR, windows = 1
#' - Moderate (default-ish): spatial_sigma = 1.5 x min(spacing), intensity_sigma = 1.0,
#'   temporal_sigma = 1 x TR, windows = 1-2
#' - Strong: spatial_sigma = 2 x min(spacing), intensity_sigma = 1.5,
#'   temporal_sigma = 1.5 x TR, windows = 2
#'
#' Tip: If your time axis has known TR, pass it via \code{temporal_spacing}.
#' For NIfTI inputs, you can get TR via:
#' \preformatted{
#'   hdr <- read_header(nifti_path)
#'   tr  <- hdr@header$pixdim[5]
#'   out <- bilateral_filter_4d(vec, mask, temporal_spacing = tr)
#' }
#'
#' @return A \code{\linkS4class{NeuroVec}} with filtered data.
#'
#' @examples
#' \donttest{
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' out  <- bilateral_filter_4d(vec, mask,
#'                             spatial_sigma = 2, intensity_sigma = 1,
#'                             temporal_sigma = 1, spatial_window = 1,
#'                             temporal_window = 1, temporal_spacing = 1)
#' }
#'
#' @seealso \code{\link{bilateral_filter}}, \code{\link{NeuroVec-class}}, \code{\link{NeuroVol-class}}
#' @export
bilateral_filter_4d <- function(vec,
                                mask,
                                spatial_sigma = 2,
                                intensity_sigma = 1,
                                temporal_sigma = 1,
                                spatial_window = 1,
                                temporal_window = 1,
                                temporal_spacing = 1,
                                range_scale = NULL) {

  if (!inherits(vec, "NeuroVec")) {
    cli::cli_abort("{.arg vec} must be a {.cls NeuroVec} object.")
  }
  if (!is.numeric(spatial_window) || spatial_window < 1) {
    cli::cli_abort("{.arg spatial_window} must be a numeric value >= 1, not {.val {spatial_window}}.")
  }
  if (!is.numeric(temporal_window) || temporal_window < 0) {
    cli::cli_abort("{.arg temporal_window} must be a numeric value >= 0, not {.val {temporal_window}}.")
  }
  if (spatial_sigma <= 0) {
    cli::cli_abort("{.arg spatial_sigma} must be positive, not {.val {spatial_sigma}}.")
  }
  if (intensity_sigma <= 0) {
    cli::cli_abort("{.arg intensity_sigma} must be positive, not {.val {intensity_sigma}}.")
  }
  if (temporal_sigma <= 0) {
    cli::cli_abort("{.arg temporal_sigma} must be positive, not {.val {temporal_sigma}}.")
  }
  if (temporal_spacing <= 0) {
    cli::cli_abort("{.arg temporal_spacing} must be positive, not {.val {temporal_spacing}}.")
  }
  if (!is.null(range_scale)) {
    if (!is.numeric(range_scale) || length(range_scale) != 1 || !is.finite(range_scale) || range_scale <= 0) {
      cli::cli_abort("{.arg range_scale} must be {.code NULL} or a single positive finite number.")
    }
  }

  # Determine mask and target space
  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vec)[1:3]))
    target_space <- space(vec)
  } else {
    if (!inherits(mask, "NeuroVol") && !inherits(mask, "LogicalNeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls NeuroVol} or {.cls LogicalNeuroVol} object.")
    }
    if (!all(dim(mask) == dim(vec)[1:3])) {
      cli::cli_abort("{.arg mask} spatial dimensions {.val {dim(mask)}} must match spatial dims of {.arg vec} {.val {dim(vec)[1:3]}}.")
    }
    if (!all(spacing(mask) == spacing(vec)[1:3])) {
      cli::cli_abort("{.arg mask} and {.arg vec} must have identical spatial spacing.")
    }
    mask.idx <- which(mask != 0)
    # Keep original 4D space
    target_space <- space(vec)
  }

  # Assemble spacing with temporal component for the 4D kernel
  sp4 <- c(spacing(vec)[1:3], temporal_spacing)

  arr <- as.array(vec)
  farr <- bilateral_filter_4d_cpp_par(arr,
                                       as.integer(mask.idx),
                                       as.integer(spatial_window),
                                       as.integer(temporal_window),
                                       spatial_sigma,
                                       intensity_sigma,
                                       temporal_sigma,
                                       sp4,
                                       if (is.null(range_scale)) NA_real_ else range_scale)

  DenseNeuroVec(farr, target_space)
}

#' Laplacian Enhancement Filter for Volumetric Images
#'
#' @description
#' This function applies a multi-layer Laplacian enhancement filter to a volumetric image (3D brain MRI data).
#' The filter enhances details while preserving edges using a non-local means approach with multiple scales.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be enhanced.
#' @param mask A \code{\linkS4class{LogicalNeuroVol}} object specifying the region to process. If not provided,
#'   the entire volume will be processed.
#' @param k An integer specifying the number of layers in the decomposition (default is 2).
#' @param patch_size An integer specifying the size of patches for non-local means. Must be odd (default is 3).
#' @param search_radius An integer specifying the radius of the search window (default is 2).
#' @param h A numeric value controlling the filtering strength. Higher values mean more smoothing (default is 0.7).
#' @param mapping_params An optional list of parameters for the enhancement mappings.
#' @param use_normalization_free Logical indicating whether to use normalization-free weights (default is TRUE).
#'
#' @return A \code{\linkS4class{NeuroVol}} object representing the enhanced image.
#'
#' @export
laplace_enhance <- function(vol, mask, k = 2, patch_size = 3, search_radius = 2,
                          h = 0.7, mapping_params = NULL,
                          use_normalization_free = TRUE) {

  if (!inherits(vol, "NeuroVol")) {
    cli::cli_abort("{.arg vol} must be a {.cls NeuroVol} object.")
  }
  if (k < 1) {
    cli::cli_abort("{.arg k} must be >= 1, not {.val {k}}.")
  }
  if (patch_size < 3 || patch_size %% 2 != 1) {
    cli::cli_abort("{.arg patch_size} must be an odd integer >= 3, not {.val {patch_size}}.")
  }
  if (search_radius < 1) {
    cli::cli_abort("{.arg search_radius} must be >= 1, not {.val {search_radius}}.")
  }
  if (h <= 0) {
    cli::cli_abort("{.arg h} must be positive, not {.val {h}}.")
  }

  # Create default mask if not provided
  if (missing(mask)) {
    mask <- LogicalNeuroVol(array(TRUE, dim(vol)), space(vol))
  } else {
    if (!inherits(mask, "LogicalNeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls LogicalNeuroVol} object.")
    }
  }

  # Call C++ implementation
  farr <- fast_multilayer_laplacian_enhancement_masked(
    as.array(vol),
    as.logical(mask),
    as.integer(k),
    as.integer(patch_size),
    as.integer(search_radius),
    h,
    mapping_params,
    use_normalization_free
  )

  # Return enhanced volume
  out <- NeuroVol(farr, space(vol))
  out
}

# ---------------------------------------------------------------------------
# Visualization enhancement for unsmoothed statistical maps
# ---------------------------------------------------------------------------

#' Selective median despike of impulsive ("salt-and-pepper") voxels
#'
#' Internal helper. Flags voxels whose deviation from the local box mean exceeds
#' \code{k} robust (MAD) scales, then replaces \emph{only} those voxels with the
#' median of their in-mask, non-flagged neighbours. Unflagged voxels keep their
#' exact value, so genuine clusters are untouched. The cheap mean-based flagging
#' is robust to isolated impulses because a single spike shifts a 27-neighbour
#' mean by only ~1/27 of its amplitude while the spike itself deviates fully.
#'
#' @param arr A numeric 3D array (carrying a \code{dim} attribute).
#' @param mask.idx Integer linear indices of in-mask voxels.
#' @param d Integer vector of array dimensions.
#' @param k Robust-deviation threshold for flagging impulses.
#' @param window Half-width (in voxels) of the neighbourhood.
#' @return A despiked copy of \code{arr}.
#' @keywords internal
#' @noRd
.despike_impulses <- function(arr, mask.idx, d, k = 3.5, window = 1L) {
  idx <- as.integer(mask.idx)
  m <- box_blur(arr, idx, as.integer(window))
  dev <- arr - m
  scale <- stats::mad(dev[mask.idx], na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    return(arr)
  }

  inmask <- logical(prod(d))
  inmask[mask.idx] <- TRUE
  flagged <- which(inmask & abs(dev) > k * scale)
  if (length(flagged) == 0L) {
    return(arr)
  }

  flagset <- logical(prod(d))
  flagset[flagged] <- TRUE

  coords <- arrayInd(flagged, .dim = d)
  offs <- as.matrix(expand.grid(seq.int(-window, window),
                                seq.int(-window, window),
                                seq.int(-window, window)))
  offs <- offs[!(offs[, 1] == 0 & offs[, 2] == 0 & offs[, 3] == 0), , drop = FALSE]

  nb_vals <- matrix(NA_real_, nrow = nrow(coords), ncol = nrow(offs))
  for (o in seq_len(nrow(offs))) {
    nc <- sweep(coords, 2L, offs[o, ], "+")
    valid <- nc[, 1] >= 1 & nc[, 1] <= d[1] &
             nc[, 2] >= 1 & nc[, 2] <= d[2] &
             nc[, 3] >= 1 & nc[, 3] <= d[3]
    if (!any(valid)) next
    lin <- nc[valid, 1] + (nc[valid, 2] - 1L) * d[1] + (nc[valid, 3] - 1L) * d[1] * d[2]
    # only borrow from in-mask voxels that are not themselves impulses
    ok <- inmask[lin] & !flagset[lin]
    rows <- which(valid)[ok]
    nb_vals[rows, o] <- arr[lin[ok]]
  }

  med <- apply(nb_vals, 1L, stats::median, na.rm = TRUE)
  has <- rowSums(!is.na(nb_vals)) > 0L
  out <- arr
  out[flagged[has]] <- med[has]
  out
}

#' Guided-filter base layer (mask-aware), returning the smoothed volume.
#'
#' @param arr Numeric 3D array with \code{dim} attribute.
#' @param idx Integer linear indices of in-mask voxels.
#' @param radius Box radius for the local linear model.
#' @param epsilon Regularisation (~ noise variance) separating noise from signal.
#' @return A numeric 3D array (the base layer).
#' @keywords internal
#' @noRd
.guided_base <- function(arr, idx, radius, epsilon) {
  mean_I  <- box_blur(arr, idx, radius)
  mean_II <- box_blur(arr * arr, idx, radius)
  var_I   <- mean_II - mean_I * mean_I
  var_I[var_I < 0] <- 0
  a <- var_I / (var_I + epsilon)
  b <- (1 - a) * mean_I
  mean_a <- box_blur(a, idx, radius)
  mean_b <- box_blur(b, idx, radius)
  mean_a * arr + mean_b
}

#' Enhance an unsmoothed statistical map for visualization
#'
#' @description
#' Removes isolated speckles from a noisy ("salt-and-pepper") statistical map so
#' it renders cleanly as an overlay, \emph{without} the peak-depressing and
#' cluster-smearing that a
#' plain Gaussian blur causes. This is a display-oriented preprocessor: it is not
#' intended to produce maps for statistical inference (it deliberately sharpens
#' for visual punch and does not preserve the null distribution).
#'
#' @details
#' The enhancement runs in three stages, each reusing the package's existing
#' filtering primitives:
#'
#' \enumerate{
#'   \item \strong{Selective despike.} Isolated impulsive voxels (flagged by a
#'     robust deviation from their local mean) are replaced by the median of
#'     their in-mask, non-impulse neighbours. All other voxels keep their exact
#'     value, so true clusters are preserved.
#'   \item \strong{Edge-preserving base smooth.} An edge-preserving filter
#'     (\code{"guided"} by default; see \code{\link{guided_filter}} /
#'     \code{\link{bilateral_filter}} / \code{\link{gaussian_blur}}) removes
#'     residual low-amplitude grain. With the guided filter, in-cluster voxels
#'     (local variance \eqn{\gg} noise variance) are returned at essentially
#'     full amplitude, so peaks are not depressed.
#'   \item \strong{Signal-gated unsharp.} Local contrast is boosted by
#'     \code{detail_gain}, weighted per-voxel by a guided-style signal weight
#'     \eqn{w = \mathrm{var}/(\mathrm{var} + \sigma^2_{noise})}. Because
#'     \eqn{w \to 0} in flat/noisy regions, the despeckled background stays
#'     smooth while real clusters gain crisp, punchy edges.
#'   }
#'
#' The noise scale \eqn{\sigma_{noise}} is estimated robustly (MAD of the
#' high-frequency residual within the mask) and is used both for the default
#' guided \code{epsilon} and for the signal weight, so the method adapts to the
#' scale of the input map (e.g. a z-map versus a raw beta map).
#'
#' @param vol A \code{\linkS4class{NeuroVol}} statistical map to enhance.
#' @param mask An optional \code{\linkS4class{NeuroVol}} / \code{\linkS4class{LogicalNeuroVol}}
#'   restricting the region processed. If omitted, all non-zero, finite voxels
#'   of \code{vol} are used.
#' @param despike Logical; run the selective median despike stage (default \code{TRUE}).
#' @param despike_k Robust-deviation threshold for flagging impulses. Lower is
#'   more aggressive (default 3.5).
#' @param despike_window Integer half-width of the despike neighbourhood (default 1L => 3x3x3).
#' @param method Edge-preserving base filter: one of \code{"guided"} (default),
#'   \code{"bilateral"}, or \code{"gaussian"}.
#' @param radius Integer box radius for the base filter and signal-weight
#'   estimation (default 2).
#' @param epsilon Guided-filter regularisation. If \code{NULL} (default) it is
#'   set to the estimated noise variance. Only used when \code{method = "guided"}.
#' @param spatial_sigma Spatial Gaussian SD for \code{method = "bilateral"} or
#'   \code{"gaussian"} (default 2).
#' @param intensity_sigma Intensity Gaussian SD for \code{method = "bilateral"} (default 1).
#' @param detail_gain Strength of the signal-gated unsharp stage. \code{0} gives
#'   a pure denoise; \code{1} restores original detail in signal regions; values
#'   \code{> 1} sharpen for display (default 1.5). Set to 0 to disable.
#' @param verbose Logical; report the number of despiked voxels and the noise
#'   estimate (default \code{FALSE}).
#'
#' @return A \code{\linkS4class{NeuroVol}} with the enhanced map (zero outside the mask).
#'
#' @examples
#' \donttest{
#' mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
#' # A noisy synthetic stat map on the mask grid
#' set.seed(1)
#' noisy <- mask * 0
#' idx <- which(mask > 0)
#' vals <- rnorm(length(idx))
#' vals[sample(length(idx), length(idx) %/% 20)] <- 8  # salt-and-pepper spikes
#' noisy[idx] <- vals
#' stat <- NeuroVol(noisy, space(mask))
#'
#' clean <- enhance_stat_map(stat, mask)
#' }
#'
#' @seealso \code{\link{guided_filter}}, \code{\link{bilateral_filter}},
#'   \code{\link{gaussian_blur}}, \code{\link{plot_overlay}}, \code{\link{plot_ortho}}
#' @export
enhance_stat_map <- function(vol, mask,
                             despike = TRUE,
                             despike_k = 3.5,
                             despike_window = 1L,
                             method = c("guided", "bilateral", "gaussian"),
                             radius = 2,
                             epsilon = NULL,
                             spatial_sigma = 2,
                             intensity_sigma = 1,
                             detail_gain = 1.5,
                             verbose = FALSE) {
  if (!inherits(vol, "NeuroVol")) {
    cli::cli_abort("{.arg vol} must be a {.cls NeuroVol} object.")
  }
  method <- match.arg(method)
  if (!is.numeric(despike_k) || length(despike_k) != 1L || despike_k <= 0) {
    cli::cli_abort("{.arg despike_k} must be a single positive number.")
  }
  if (!is.numeric(radius) || length(radius) != 1L || radius < 1) {
    cli::cli_abort("{.arg radius} must be a single number >= 1, not {.val {radius}}.")
  }
  if (!is.numeric(detail_gain) || length(detail_gain) != 1L || detail_gain < 0) {
    cli::cli_abort("{.arg detail_gain} must be a single non-negative number.")
  }
  for (nm in c("spatial_sigma", "intensity_sigma")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      cli::cli_abort("{.arg {nm}} must be a single finite positive number, not {.val {v}}.")
    }
  }
  if (!is.null(epsilon) &&
      (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0)) {
    cli::cli_abort("{.arg epsilon} must be {.code NULL} or a single positive finite number.")
  }

  d <- dim(vol)
  arr <- as.array(vol)
  arr[!is.finite(arr)] <- 0

  if (missing(mask)) {
    mask.idx <- which(arr != 0)
    target_space <- space(vol)
  } else {
    if (!inherits(mask, "NeuroVol")) {
      cli::cli_abort("{.arg mask} must be a {.cls NeuroVol} object.")
    }
    if (!identical(as.integer(dim(mask)), as.integer(d))) {
      cli::cli_abort(c(
        "{.arg mask} and {.arg vol} must have the same dimensions.",
        i = "{.arg vol} is {.val {d}}, {.arg mask} is {.val {dim(mask)}}."
      ))
    }
    mask.idx <- which(mask != 0)
    target_space <- space(mask)
  }
  if (length(mask.idx) == 0L) {
    cli::cli_abort("Mask is empty; nothing to enhance.")
  }
  idx <- as.integer(mask.idx)
  radius <- as.integer(round(radius))

  # ---- Stage 1: selective median despike ------------------------------------
  despiked <- arr
  n_despiked <- 0L
  if (isTRUE(despike)) {
    despiked <- .despike_impulses(arr, mask.idx, d, k = despike_k,
                                  window = as.integer(despike_window))
    n_despiked <- sum(despiked[mask.idx] != arr[mask.idx])
  }

  # ---- Robust noise estimate (MAD of high-frequency residual) ---------------
  hf <- despiked - box_blur(despiked, idx, 1L)
  noise_sd <- stats::mad(hf[mask.idx], na.rm = TRUE)
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    noise_sd <- stats::sd(hf[mask.idx], na.rm = TRUE)
  }
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    noise_sd <- 1
  }
  noise_var <- noise_sd^2

  # ---- Stage 2: edge-preserving base layer ----------------------------------
  base <- switch(method,
    guided = .guided_base(despiked, idx, radius,
                          if (is.null(epsilon)) noise_var else epsilon),
    bilateral = bilateral_filter_cpp(despiked, idx, max(1L, radius),
                                     spatial_sigma, intensity_sigma,
                                     spacing(vol), NA_real_),
    gaussian = .gaussian_blur_engine(despiked, idx, max(1L, radius),
                                     spatial_sigma, spacing(vol))
  )

  # ---- Stage 3: signal-gated unsharp ----------------------------------------
  out <- base
  if (detail_gain > 0) {
    mean_I <- box_blur(despiked, idx, radius)
    mean_II <- box_blur(despiked * despiked, idx, radius)
    var_I <- mean_II - mean_I * mean_I
    var_I[var_I < 0] <- 0
    w <- var_I / (var_I + noise_var)          # guided-style signal weight in [0,1]

    radius_broad <- as.integer(max(radius + 1L, ceiling(radius * 2)))
    broad <- box_blur(base, idx, radius_broad)
    out <- base + detail_gain * w * (base - broad)
  }

  if (isTRUE(verbose)) {
    cli::cli_inform(c(
      "i" = "enhance_stat_map: {n_despiked} voxel{?s} despiked.",
      "i" = "Noise scale (MAD) = {.val {round(noise_sd, 4)}}; method = {.val {method}}."
    ))
  }

  res <- numeric(prod(d))
  res[mask.idx] <- out[mask.idx]
  dim(res) <- d
  NeuroVol(res, target_space)
}
