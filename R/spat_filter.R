

#' Blur a volumetric image with an isotropic discrete Gaussian kernel
#'
#' This function smooths a volumetric image (3D brain MRI data) by applying an isotropic discrete Gaussian kernel.
#' The blurring is performed within the specified image mask, and the kernel's standard deviation and window size can be customized.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.
#' @param mask An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask that defines the region where the blurring is applied. If not provided, the entire volume is considered.
#' @param sigma A numeric value specifying the standard deviation of the Gaussian kernel (default is 2).
#' @param window An integer specifying the number of voxels around the center voxel to include on each side. For example, window=1 for a 3x3x3 kernel (default is 1).
#'
#' @return A smoothed image of class \code{\linkS4class{NeuroVol}}.
#'
#' @examples
#' brain_mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#'
#' # Apply Gaussian blurring to the brain volume
#' blurred_vol <- gaussian_blur(brain_mask, brain_mask, sigma = 2, window = 1)
#'
#' @export
gaussian_blur <- function(vol, mask, sigma=2, window=1) {
  assert_that(window >= 1)
  assert_that(sigma > 0)

  if (missing(mask)) {
    mask.idx <- 1:prod(dim(vol))
  } else {
    mask.idx <- which(mask!=0)
  }

  arr <- as.array(vol)
  farr <- gaussian_blur_cpp(arr, as.integer(mask.idx), as.integer(window), sigma, spacing(vol))

  out <- NeuroVol(farr, space(mask))
  out
}

#' Filter a volumetric image with an edge-preserving "guided" filter
#'
#' This function applies a guided filter to a volumetric image (3D brain MRI data) to perform edge-preserving smoothing.
#' The guided filter is an edge-preserving filter that smooths the image while preserving the edges, providing a balance between noise reduction and edge preservation.
#'
#' @param vol A \code{\linkS4class{NeuroVol}} object representing the image volume to be filtered.
#' @param radius An integer specifying the spatial radius of the filter (default is 4).
#' @param epsilon A numeric value specifying the variance constant, which controls the degree of smoothing (default is .7^2).
#'
#' @return A filtered image of class \code{\linkS4class{NeuroVol}}.
#'
#' @examples
#' # Load an example brain volume
#' brain_vol <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
#'
#' # Apply guided filtering to the brain volume
#' \dontrun{
#' filtered_vol <- guided_filter(brain_vol, radius = 4, epsilon = .7^2)
#' }
#'
#' @references
#' Guided Image Filtering: Kaiming He, Jian Sun, and Xiaoou Tang, "Guided Image Filtering," IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 35, No. 6, pp. 1397-1409, June, 2013.
#' https://en.wikipedia.org/wiki/Guided_filter
#'
#' @export
guided_filter <- function(vol, radius=4, epsilon=.7^2) {
  # pset <- patch_set(vol, c(3,3,3))
  #
  # ab <- do.call(rbind, parallel::mclapply(pset, function(v) {
  #   mu_k <- mean(v)
  #   ak <- sum(((v^2 - mu_k^2)))/length(v)
  #   ak <- ak/(sd(v) + epsilon)
  #   bk <- mu_k - ak*mu_k
  #   c(ak, bk)
  # }, mc.cores=ncores))
  #
  # resd <- do.call(rbind, parallel::mclapply(1:length(pset), function(i) {
  #   v <- pset[[i]]
  #   qi <- ab[i,1] * v + ab[i,2]
  #   data.frame(q=qi, idx=attr(v, "idx"))
  # }, mc.cores=ncores))
  #
  # resq <- resd %>% group_by(idx) %>% summarize(q=mean(q))
  # out <- NeuroVol(resq$q, space(vol), indices=resq$idx)
  # out

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
