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
#' brain_mask <- read_vol(system.file("extdata", "global_mask.nii", package = "neuroim2"))
#'
#' # Apply Gaussian blurring to the brain volume
#' blurred_vol <- gaussian_blur(brain_mask, brain_mask, sigma = 2, window = 1)
#'
#' # View a slice of the original and blurred volumes
#' image(brain_mask[,,50])
#' image(blurred_vol[,,50])
#'
#' @seealso 
#' \code{\link{NeuroVol-class}}, \code{\link{LogicalNeuroVol-class}}, \code{\link{bilateral_filter}}
#'
#' @references
#' Gaussian blur: https://en.wikipedia.org/wiki/Gaussian_blur
#'
#' @export
gaussian_blur <- function(vol, mask, sigma = 2, window = 1) {
  assert_that(window >= 1)
  assert_that(sigma > 0)

  if (missing(mask)) {
    mask.idx <- 1:prod(dim(vol))
  } else {
    mask.idx <- which(mask != 0)
  }

  arr <- as.array(vol)
  farr <- gaussian_blur_cpp(arr, as.integer(mask.idx), as.integer(window), sigma, spacing(vol))

  out <- NeuroVol(farr, space(mask))
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
#' brain_vol <- read_vol(system.file("extdata", "global_mask.nii", package = "neuroim2"))
#'
#' # Apply guided filtering to the brain volume
#' \dontrun{
#' filtered_vol <- guided_filter(brain_vol, radius = 4, epsilon = 0.49)
#'
#' # Visualize a slice of the original and filtered volumes
#' par(mfrow = c(1, 2))
#' image(brain_vol[,,50], main = "Original")
#' image(filtered_vol[,,50], main = "Filtered")
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
#' brain_mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
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

  if (missing(mask)) {
    mask.idx <- 1:prod(dim(vol))
  } else {
    mask.idx <- which(mask!=0)
  }

  arr <- as.array(vol)
  farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window), spatial_sigma, intensity_sigma, spacing(vol))

  out <- NeuroVol(farr, space(mask))
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
#' brain_mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#' vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' out <- bilateral_filter_vec(vec,brain_mask)
#' @noRd
bilateral_filter_vec <- function(vec, mask, spatial_sigma=2, intensity_sigma=1, window=1) {
  assert_that(inherits(vec, "NeuroVec"))
  assert_that(window >= 1)
  assert_that(spatial_sigma > 0)
  assert_that(intensity_sigma > 0)

  if (missing(mask)) {
    mask.idx <- 1:prod(dim(mask))
  } else {
    mask.idx <- which(mask!=0)
  }

  res<- lapply(1:dim(vec)[4], function(i) {
    #print(i)
    arr <- as.array(vec[[i]])
    farr <- bilateral_filter_cpp(arr, as.integer(mask.idx), as.integer(window), spatial_sigma, intensity_sigma, spacing(vec)[1:3])
    NeuroVol(farr, space(mask))
  })

  do.call(concat,res)

}


