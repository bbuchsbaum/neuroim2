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
#' function for efficiency. Set any FWHM argument to \code{0} to disable that
#' smoothing step.
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
  if (!inherits(mask, "NeuroVol")) {
    cli::cli_abort("{.arg mask} must inherit from {.cls NeuroVol}.")
  }

  assert_scalar <- function(x, name, whole = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
      cli::cli_abort("{.arg {name}} must be a single finite numeric value.")
    }
    if (whole && x != floor(x)) {
      cli::cli_abort("{.arg {name}} must be a whole number.")
    }
    invisible(x)
  }

  assert_scalar(n_time, "n_time", whole = TRUE)
  assert_scalar(TR, "TR")
  assert_scalar(spatial_fwhm, "spatial_fwhm")
  assert_scalar(ar_mean, "ar_mean")
  assert_scalar(ar_sd, "ar_sd")
  assert_scalar(noise_sd, "noise_sd")
  assert_scalar(hetero_fwhm, "hetero_fwhm")
  assert_scalar(hetero_strength, "hetero_strength")
  assert_scalar(global_amp, "global_amp")
  assert_scalar(global_rho, "global_rho")
  assert_scalar(n_factors, "n_factors", whole = TRUE)
  assert_scalar(factor_fwhm, "factor_fwhm")
  assert_scalar(factor_rho, "factor_rho")

  if (n_time < 1L) cli::cli_abort("{.arg n_time} must be >= 1.")
  if (TR <= 0) cli::cli_abort("{.arg TR} must be > 0.")
  if (spatial_fwhm < 0) cli::cli_abort("{.arg spatial_fwhm} must be >= 0.")
  if (hetero_fwhm < 0) cli::cli_abort("{.arg hetero_fwhm} must be >= 0.")
  if (factor_fwhm < 0) cli::cli_abort("{.arg factor_fwhm} must be >= 0.")
  if (ar_mean < 0 || ar_mean >= 1) cli::cli_abort("{.arg ar_mean} must be in [0, 1).")
  if (ar_sd < 0) cli::cli_abort("{.arg ar_sd} must be >= 0.")
  if (noise_sd <= 0) cli::cli_abort("{.arg noise_sd} must be > 0.")
  if (global_amp < 0) cli::cli_abort("{.arg global_amp} must be >= 0.")
  if (global_rho < 0 || global_rho >= 1) cli::cli_abort("{.arg global_rho} must be in [0, 1).")
  if (n_factors < 0L) cli::cli_abort("{.arg n_factors} must be >= 0.")
  if (factor_rho < 0 || factor_rho >= 1) cli::cli_abort("{.arg factor_rho} must be in [0, 1).")
  if (!is.logical(return_centered) || length(return_centered) != 1L || is.na(return_centered)) {
    cli::cli_abort("{.arg return_centered} must be TRUE or FALSE.")
  }
  
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
  spatial_window <- if (spatial_sigma > 0) max(1L, ceiling(3 * spatial_sigma / min(voxsize))) else 0L
  hetero_window <- if (hetero_sigma > 0) max(1L, ceiling(3 * hetero_sigma / min(voxsize))) else 0L
  factor_window <- if (factor_sigma > 0) max(1L, ceiling(3 * factor_sigma / min(voxsize))) else 0L
  
  # Create mask as LogicalNeuroVol for gaussian_blur
  mask_vol <- as.logical(mask)

  smooth_vol <- function(vol, sigma, window) {
    if (sigma > 0) {
      gaussian_blur(vol, mask_vol, sigma = sigma, window = as.integer(window))
    } else {
      vol
    }
  }
  
  # Pre-compute scaling factor for spatial smoothing
  # Create calibration noise and smooth it to determine scaling
  calib_array <- array(0, dm)
  calib_array[mask_idx] <- rnorm(n_vox)
  calib_vol <- NeuroVol(calib_array, space(mask))
  calib_smooth <- smooth_vol(calib_vol, spatial_sigma, spatial_window)
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
                                window = max(1L, as.integer(ceiling(spatial_window / 4))))
    rho_array <- as.array(rho_smooth)
  }
  
  rho_vec <- pmax(pmin(rho_array[mask_idx], 0.95), 0.0)
  
  # Generate heteroscedastic spatial SD field
  het_array <- array(0, dm)
  het_array[mask_idx] <- rnorm(n_vox)
  het_vol <- NeuroVol(het_array, space(mask))
  het_smooth <- smooth_vol(het_vol, hetero_sigma, hetero_window)
  het_field <- as.array(het_smooth)[mask_idx]
  if (hetero_strength == 0) {
    het_field <- numeric(n_vox)
  } else {
    het_mean <- mean(het_field)
    het_sd <- sd(het_field)
    if (is.finite(het_sd) && het_sd > 0) {
      het_field <- as.numeric((het_field - het_mean) / het_sd)
    } else {
      het_field <- numeric(n_vox)
    }
  }
  sd_vox <- as.numeric(noise_sd * exp(hetero_strength * het_field))
  
  # AR(1) innovation SD per voxel
  sigma_e <- sd_vox * sqrt(pmax(1e-8, 1 - rho_vec^2))
  
  # Generate optional global signal
  g <- numeric(n_time)
  if (global_amp > 0) {
    g_e_sd <- sqrt(1 - global_rho^2)
    g[1] <- rnorm(1, sd = 1)
    if (n_time >= 2L) {
      for (t in 2:n_time) {
        g[t] <- global_rho * g[t-1] + rnorm(1, sd = g_e_sd)
      }
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
      m_smooth <- smooth_vol(m_vol, factor_sigma, factor_window)
      v <- as.array(m_smooth)[mask_idx]
      v <- v / sqrt(sum(v^2) + 1e-12)  # L2 normalize
      F_maps[, k] <- v
    }
    
    # Generate AR(1) time courses for factors
    Z_tc <- matrix(0, nrow = n_factors, ncol = n_time)
    e_sd <- sqrt(1 - factor_rho^2)
    Z_tc[, 1] <- rnorm(n_factors, sd = 1)
    if (n_time >= 2L) {
      for (t in 2:n_time) {
        Z_tc[, t] <- factor_rho * Z_tc[, t-1] + rnorm(n_factors, sd = e_sd)
      }
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
    eps_smooth <- smooth_vol(eps_vol, spatial_sigma, spatial_window)
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
    out_mat <- matrix(out, nrow = prod(dm), ncol = n_time)
    out_mat[mask_idx, ] <- out_mat[mask_idx, , drop = FALSE] -
      rowMeans(out_mat[mask_idx, , drop = FALSE])
    out <- array(out_mat, dim = c(dm, n_time))
  }
  
  # Create 4D NeuroSpace by adding time dimension to mask space
  space_4d <- add_dim(space(mask), n_time)
  
  # Return as NeuroVec
  result <- NeuroVec(out, space_4d)
  attr(result, "TR") <- TR
  result
}
