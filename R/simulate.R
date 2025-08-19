#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Simulate fMRI Data
#'
#' @description
#' Generates synthetic 4D fMRI data with realistic spatiotemporal properties including
#' temporal autocorrelation, spatial smoothness, heteroscedasticity, and optional 
#' global signal fluctuations and latent components.
#'
#' @param mask A \code{\linkS4class{NeuroVol}} object defining the brain mask region.
#'   Can be binary or continuous (non-zero values define the mask).
#' @param n_time Integer specifying the number of time points to simulate.
#' @param TR Numeric value for the repetition time in seconds (default = 2.0).
#'   Currently used only for metadata.
#' @param spatial_fwhm Numeric value specifying the spatial smoothness in mm 
#'   (full width at half maximum) applied to each timepoint (default = 6).
#' @param ar_mean Numeric value for the mean of the AR(1) coefficient distribution 
#'   across voxels (default = 0.45).
#' @param ar_sd Numeric value for the standard deviation of the AR(1) coefficient 
#'   distribution (default = 0.08).
#' @param noise_sd Numeric value for the baseline noise standard deviation 
#'   (default = 1.0).
#' @param hetero_fwhm Numeric value for the spatial scale (FWHM in mm) of the 
#'   heteroscedasticity field (default = 20).
#' @param hetero_strength Numeric value controlling the strength of spatial 
#'   heteroscedasticity on log scale (default = 0.6).
#' @param global_amp Numeric value for the amplitude of global signal fluctuations 
#'   as a fraction of median noise (default = 0.2). Set to 0 to disable.
#' @param global_rho Numeric value for the AR(1) coefficient of global signal 
#'   (default = 0.85).
#' @param n_factors Integer specifying the number of latent spatial components 
#'   (default = 4). Set to 0 to disable.
#' @param factor_fwhm Numeric value for the spatial smoothness (FWHM in mm) of 
#'   latent component maps (default = 12).
#' @param factor_rho Numeric value for the AR(1) coefficient of latent component 
#'   time courses (default = 0.8).
#' @param seed Integer seed for random number generation (default = NULL for no seed).
#' @param return_centered Logical indicating whether to center each voxel's time series 
#'   to mean zero (default = TRUE).
#'
#' @return A \code{\linkS4class{NeuroVec}} object containing the simulated 4D fMRI data.
#'
#' @details
#' The simulation combines several realistic features:
#' \itemize{
#'   \item Voxel-wise AR(1) temporal autocorrelation with spatial variation
#'   \item Spatial smoothing applied to innovations for realistic spatial correlation
#'   \item Heteroscedastic noise with smooth spatial modulation
#'   \item Optional low-frequency global signal fluctuations
#'   \item Optional latent spatial components resembling resting-state networks
#' }
#'
#' The spatial smoothing uses the package's optimized \code{\link{gaussian_blur}} 
#' function for efficiency.
#'
#' @examples
#' # Create a simple spherical mask
#' dims <- c(32, 32, 20)
#' mask_array <- array(FALSE, dims)
#' center <- dims / 2
#' for (i in 1:dims[1]) {
#'   for (j in 1:dims[2]) {
#'     for (k in 1:dims[3]) {
#'       if (sum(((c(i,j,k) - center) / (dims/3))^2) <= 1) {
#'         mask_array[i,j,k] <- TRUE
#'       }
#'     }
#'   }
#' }
#' 
#' mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))
#' 
#' # Simulate 100 time points
#' sim_data <- simulate_fmri(mask, n_time = 100, seed = 42)
#' 
#' # Check dimensions
#' dim(sim_data)  # Should be c(32, 32, 20, 100)
#'
#' @references
#' Welvaert, M., & Rosseel, Y. (2013). On the definition of signal-to-noise ratio 
#' and contrast-to-noise ratio for fMRI data. PloS one, 8(11), e77089.
#'
#' @export
#' @importFrom stats rnorm sd filter
simulate_fmri <- function(mask,
                         n_time,
                         TR = 2.0,
                         spatial_fwhm = 6,
                         ar_mean = 0.45,
                         ar_sd = 0.08,
                         noise_sd = 1.0,
                         hetero_fwhm = 20,
                         hetero_strength = 0.6,
                         global_amp = 0.2,
                         global_rho = 0.85,
                         n_factors = 4,
                         factor_fwhm = 12,
                         factor_rho = 0.8,
                         seed = NULL,
                         return_centered = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  stopifnot(inherits(mask, "NeuroVol"))
  stopifnot(n_time > 0)
  stopifnot(spatial_fwhm >= 0)
  stopifnot(ar_mean >= 0 && ar_mean < 1)
  stopifnot(ar_sd >= 0)
  stopifnot(noise_sd > 0)
  
  # Get mask dimensions and indices
  dm <- dim(mask)
  mask_array <- as.array(mask) != 0
  mask_idx <- which(mask_array)
  n_vox <- length(mask_idx)
  
  if (n_vox == 0) stop("Mask is empty")
  
  # Get voxel sizes from mask space
  voxsize <- spacing(mask)
  
  # Helper function to convert FWHM to sigma for gaussian_blur
  fwhm_to_sigma <- function(fwhm_mm) {
    fwhm_mm / (2 * sqrt(2 * log(2)))
  }
  
  # Convert FWHM values to sigma for gaussian_blur
  spatial_sigma <- fwhm_to_sigma(spatial_fwhm)
  hetero_sigma <- fwhm_to_sigma(hetero_fwhm)
  factor_sigma <- fwhm_to_sigma(factor_fwhm)
  
  # Compute window size for gaussian_blur (3 sigma on each side)
  spatial_window <- max(1, ceiling(3 * spatial_sigma / min(voxsize)))
  hetero_window <- max(1, ceiling(3 * hetero_sigma / min(voxsize)))
  factor_window <- max(1, ceiling(3 * factor_sigma / min(voxsize)))
  
  # Create mask as LogicalNeuroVol for gaussian_blur
  mask_vol <- as.logical(mask)
  
  # Pre-compute scaling factor for spatial smoothing
  # Create calibration noise and smooth it to determine scaling
  calib_array <- array(0, dm)
  calib_array[mask_idx] <- rnorm(n_vox)
  calib_vol <- NeuroVol(calib_array, space(mask))
  calib_smooth <- gaussian_blur(calib_vol, mask_vol, sigma = spatial_sigma, window = spatial_window)
  sm_sd <- sd(as.array(calib_smooth)[mask_idx])
  smooth_scale <- if (!is.na(sm_sd) && sm_sd > 0) 1 / sm_sd else 1
  
  # Generate voxel-wise AR(1) parameters
  raw_rho <- pmax(pmin(rnorm(n_vox, mean = ar_mean, sd = ar_sd), 0.95), 0.0)
  rho_array <- array(0, dm)
  rho_array[mask_idx] <- raw_rho
  
  # Optionally smooth the rho field slightly for spatial coherence
  if (spatial_sigma > 0) {
    rho_vol <- NeuroVol(rho_array, space(mask))
    rho_smooth <- gaussian_blur(rho_vol, mask_vol, 
                                sigma = spatial_sigma / 4,  # Light smoothing
                                window = max(1, spatial_window / 4))
    rho_array <- as.array(rho_smooth)
  }
  
  rho_vec <- pmax(pmin(rho_array[mask_idx], 0.95), 0.0)
  
  # Generate heteroscedastic spatial SD field
  het_array <- array(0, dm)
  het_array[mask_idx] <- rnorm(n_vox)
  het_vol <- NeuroVol(het_array, space(mask))
  het_smooth <- gaussian_blur(het_vol, mask_vol, sigma = hetero_sigma, window = hetero_window)
  het_field <- as.array(het_smooth)[mask_idx]
  het_field <- as.numeric(base::scale(het_field))  # Standardize
  sd_vox <- as.numeric(noise_sd * exp(hetero_strength * het_field))
  
  # AR(1) innovation SD per voxel
  sigma_e <- sd_vox * sqrt(pmax(1e-8, 1 - rho_vec^2))
  
  # Generate optional global signal
  g <- numeric(n_time)
  if (global_amp > 0) {
    g_e_sd <- sqrt(1 - global_rho^2)
    g[1] <- rnorm(1, sd = 1)
    for (t in 2:n_time) {
      g[t] <- global_rho * g[t-1] + rnorm(1, sd = g_e_sd)
    }
    g <- g * (global_amp * stats::median(sd_vox))
  }
  
  # Generate optional latent factor components
  F_maps <- NULL
  Z_tc <- NULL
  if (n_factors > 0) {
    F_maps <- matrix(0, nrow = n_vox, ncol = n_factors)
    for (k in seq_len(n_factors)) {
      m_array <- array(0, dm)
      m_array[mask_idx] <- rnorm(n_vox)
      m_vol <- NeuroVol(m_array, space(mask))
      m_smooth <- gaussian_blur(m_vol, mask_vol, sigma = factor_sigma, window = factor_window)
      v <- as.array(m_smooth)[mask_idx]
      v <- v / sqrt(sum(v^2) + 1e-12)  # L2 normalize
      F_maps[, k] <- v
    }
    
    # Generate AR(1) time courses for factors
    Z_tc <- matrix(0, nrow = n_factors, ncol = n_time)
    e_sd <- sqrt(1 - factor_rho^2)
    Z_tc[, 1] <- rnorm(n_factors, sd = 1)
    for (t in 2:n_time) {
      Z_tc[, t] <- factor_rho * Z_tc[, t-1] + rnorm(n_factors, sd = e_sd)
    }
    Z_tc <- Z_tc * (0.3 * stats::median(sd_vox))
  }
  
  # Initialize AR(1) state from stationary distribution
  state <- rnorm(n_vox, sd = sd_vox)
  
  # Create output 4D array
  out <- array(0, dim = c(dm, n_time))
  
  # Main simulation loop
  for (t in seq_len(n_time)) {
    # Generate spatially correlated innovation
    eps_array <- array(0, dm)
    eps_array[mask_idx] <- rnorm(n_vox)
    eps_vol <- NeuroVol(eps_array, space(mask))
    eps_smooth <- gaussian_blur(eps_vol, mask_vol, sigma = spatial_sigma, window = spatial_window)
    e_vec <- as.array(eps_smooth)[mask_idx] * smooth_scale * sigma_e
    
    # AR(1) update
    state <- rho_vec * state + e_vec
    
    # Add global signal
    y <- state
    if (global_amp > 0) {
      y <- y + g[t]
    }
    
    # Add latent factors
    if (n_factors > 0) {
      y <- y + as.numeric(F_maps %*% Z_tc[, t])
    }
    
    # Store in output array
    vol_t <- array(0, dm)
    vol_t[mask_idx] <- y
    out[,,,t] <- vol_t
  }
  
  # Optionally center each voxel time series
  if (return_centered) {
    mu <- apply(out, 1:3, mean)
    for (t in seq_len(n_time)) {
      out[,,,t] <- out[,,,t] - mu
    }
  }
  
  # Create 4D NeuroSpace by adding time dimension to mask space
  space_4d <- add_dim(space(mask), n_time)
  
  # Return as NeuroVec
  result <- NeuroVec(out, space_4d)
  attr(result, "TR") <- TR
  result
}