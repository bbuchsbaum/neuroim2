

#' Blur a volumetric image with an isotropic discrete Gaussian kernel
#'
#' @param vol the image volume as a \code{NeuroVol}
#' @param mask an image mask as a \code{LogicalNeuroVol}
#' @param sigma the standard deviation of the Gaussian
#' @param window the number of voxels around the center voxel to include on each side (window=1 for a 3x3x3 kernel).
#'
#' @return a smoothed image of class \code{NeuroVol}
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
#' @param vol the image volume as a \code{NeuroVol}
#' @param radius the spatial radius of the filter
#' @param epsilon variance constant
#' @return a filtered image of class \code{NeuroVol}
#' @export
#' @references
#' https://en.wikipedia.org/wiki/Guided_filter
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
