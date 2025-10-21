#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Spatial Filtering Methods for Neuroimaging Data
#'
#' @name spatial-filter
#' @description Methods for applying spatial filters to neuroimaging data
NULL

#' @importFrom assertthat assert_that
#' @importFrom methods new
#' @importFrom stats dnorm
NULL

#' Gaussian Blur for Volumetric Images
#'
#' @description
#' This function applies an isotropic discrete Gaussian kernel to smooth a volumetric image (3D brain MRI data).
#' The blurring is performed within a specified image mask, with customizable kernel parameters.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask.
#'   This mask defines the region where the blurring is applied. If not provided, the entire volume is processed.
#' @param sigma A numeric value specifying the standard deviation of the Gaussian kernel. Default is 2.
#' @param window An integer specifying the kernel size. It represents the number of voxels to include
#'   on each side of the center voxel. For example, window=1 results in a 3x3x3 kernel. Default is 1.
#'
#' @return A \code{\linkS4class{NeuroVol}} object representing the smoothed image.
#'
#' @details
#' The function uses a C++ implementation for efficient Gaussian blurring. The blurring is applied
#' only to voxels within the specified mask (or the entire volume if no mask is provided).
#' The kernel size is determined by the 'window' parameter, and its shape by the 'sigma' parameter.
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
#' @seealso
#' \code{\link{NeuroVol-class}}, \code{\link{LogicalNeuroVol-class}}, \code{\link{bilateral_filter}}
#'
#' @references
#' Gaussian blur: https://en.wikipedia.org/wiki/Gaussian_blur
#'
#' @export
gaussian_blur <- function(vol, mask, sigma = 2, window = 1) {
  assert_that(inherits(vol, "NeuroVol"),
              msg = "vol must be a NeuroVol object")
  assert_that(window >= 1,
              msg = "window must be >= 1")
  assert_that(sigma > 0,
              msg = "sigma must be positive")
  if (!missing(mask)) {
    assert_that(inherits(mask, "NeuroVol"),
                msg = "mask must be a NeuroVol object")
  }

  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vol)))
    target_space <- space(vol)
  } else {
    mask.idx <- which(mask != 0)
    target_space <- space(mask)
  }

  arr <- as.array(vol)
  farr <- gaussian_blur_cpp(arr, as.integer(mask.idx), as.integer(window), sigma, spacing(vol))

  out <- NeuroVol(farr, target_space)
  out
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
  assert_that(inherits(vol, "NeuroVol"),
              msg = "vol must be a NeuroVol object")
  assert_that(radius >= 1,
              msg = "radius must be >= 1")
  assert_that(epsilon > 0,
              msg = "epsilon must be positive")

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
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask that defines the region where the filtering is applied. If not provided, the entire volume is considered.
#' @param spatial_sigma A numeric value specifying the standard deviation of the spatial Gaussian kernel (default is 2).
#' @param intensity_sigma A numeric value specifying the standard deviation of the intensity Gaussian kernel (default is 25).
#' @param window An integer specifying the number of voxels around the center voxel to include on each side. For example, window=1 for a 3x3x3 kernel (default is 1).
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
bilateral_filter <- function(vol, mask, spatial_sigma=2, intensity_sigma=1, window=1) {
  assert_that(window >= 1)
  assert_that(spatial_sigma > 0)
  assert_that(intensity_sigma > 0)
  if (!missing(mask)) {
    assert_that(inherits(mask, "NeuroVol"),
                msg = "mask must be a NeuroVol object")
  }

  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vol)))
    target_space <- space(vol)
  } else {
    mask.idx <- which(mask!=0)
    target_space <- space(mask)
  }

  arr <- as.array(vol)
  farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window), spatial_sigma, intensity_sigma, spacing(vol))

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
#' @return A NeuroVec object with the filtered volumes.
#' @examples
#' brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' out <- bilateral_filter_vec(vec,brain_mask)
#' @noRd
bilateral_filter_vec <- function(vec, mask, spatial_sigma=2, intensity_sigma=1, window=1) {
  assert_that(inherits(vec, "NeuroVec"))
  assert_that(window >= 1)
  assert_that(spatial_sigma > 0)
  assert_that(intensity_sigma > 0)

  if (missing(mask)) {
    mask.idx <- seq_len(prod(dim(vec)[1:3]))
    target_space <- space(vec[[1]])
  } else {
    assert_that(inherits(mask, "NeuroVol"),
                msg = "mask must be a NeuroVol object")
    mask.idx <- which(mask!=0)
    target_space <- space(mask)
  }

  res<- lapply(seq_len(dim(vec)[4]), function(i) {
    vol_i <- vec[[i]]
    arr <- as.array(vol_i)
    farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window), spatial_sigma, intensity_sigma, spacing(vec)[1:3])
    NeuroVol(farr, target_space)
  })

  do.call(concat,res)

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

  assert_that(inherits(vol, "NeuroVol"))
  assert_that(k >= 1)
  assert_that(patch_size >= 3 && patch_size %% 2 == 1)
  assert_that(search_radius >= 1)
  assert_that(h > 0)

  # Create default mask if not provided
  if (missing(mask)) {
    mask <- LogicalNeuroVol(array(TRUE, dim(vol)), space(vol))
  } else {
    assert_that(inherits(mask, "LogicalNeuroVol"))
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
