#' @include all_class.R
NULL

#' Downsampling Methods for Neuroimaging Objects
#'
#' @name neuro-downsample
#' @description Methods for downsampling neuroimaging objects to lower spatial resolution
NULL

#' Calculate new dimensions for downsampling
#' @keywords internal
#' @noRd
calculate_downsample_dims <- function(current_dims, current_spacing, 
                                     spacing = NULL, factor = NULL, outdim = NULL) {
  
  # current_dims should be 4D (x, y, z, t)
  spatial_dims <- current_dims[1:3]
  time_dim <- current_dims[4]
  
  if (!is.null(spacing)) {
    # Validate spacing values
    if (length(spacing) != 3) {
      stop("spacing must be a numeric vector of length 3")
    }
    if (any(spacing <= 0)) {
      stop("spacing values must be positive")
    }
    
    # Downsample to achieve target spacing
    # New dims = old dims * (old spacing / new spacing)
    new_spatial_dims <- round(spatial_dims * (current_spacing[1:3] / spacing))
    new_spatial_dims <- pmax(new_spatial_dims, 1) # Ensure at least 1 voxel
    
  } else if (!is.null(factor)) {
    # Downsample by factor
    if (length(factor) == 1) {
      factor <- rep(factor, 3)
    }
    if (any(factor <= 0) || any(factor > 1)) {
      stop("factor must be between 0 (exclusive) and 1 (inclusive)")
    }
    new_spatial_dims <- round(spatial_dims * factor)
    new_spatial_dims <- pmax(new_spatial_dims, 1)
    
  } else if (!is.null(outdim)) {
    # Check aspect ratio preservation
    if (length(outdim) != 3) {
      stop("outdim must have exactly 3 values for spatial dimensions")
    }
    
    # Calculate the implied factors
    factors <- outdim / spatial_dims
    
    # Check if aspect ratios are preserved (within tolerance)
    factor_range <- range(factors)
    mean_factor <- mean(factors)
    
    # Check for zero mean (which would happen if outdim has zeros)
    if (mean_factor == 0) {
      stop("Output dimensions cannot be zero")
    }
    
    if ((factor_range[2] - factor_range[1]) / mean_factor > 0.01) {
      warning("Output dimensions do not preserve aspect ratios. Using uniform scaling based on smallest dimension.")
      # Use the smallest factor to preserve aspect ratio
      uniform_factor <- min(factors)
      new_spatial_dims <- round(spatial_dims * uniform_factor)
    } else {
      new_spatial_dims <- outdim
    }
    new_spatial_dims <- pmax(new_spatial_dims, 1)
    
  } else {
    stop("Must specify one of: spacing, factor, or outdim")
  }
  
  # Return 4D dimensions with preserved time
  c(new_spatial_dims, time_dim)
}

#' Downsample a DenseNeuroVec
#'
#' @rdname downsample-methods
#' @param x A DenseNeuroVec object to downsample
#' @param spacing Target voxel spacing (numeric vector of length 3)
#' @param factor Downsampling factor (single value or vector of length 3, between 0 and 1)
#' @param outdim Target output dimensions (numeric vector of length 3)
#' @param method Downsampling method (currently only "box" for box averaging)
#'
#' @examples
#' # Create a sample 4D image
#' data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
#' space <- NeuroSpace(dim = c(64, 64, 32, 10), 
#'                     origin = c(0, 0, 0),
#'                     spacing = c(2, 2, 2))
#' nvec <- DenseNeuroVec(data, space)
#'
#' # Downsample by factor
#' nvec_down1 <- downsample(nvec, factor = 0.5)
#' 
#' # Downsample to target spacing
#' nvec_down2 <- downsample(nvec, spacing = c(4, 4, 4))
#' 
#' # Downsample to target dimensions
#' nvec_down3 <- downsample(nvec, outdim = c(32, 32, 16))
#'
#' @export
setMethod(f="downsample", signature=signature("DenseNeuroVec"),
          def=function(x, spacing=NULL, factor=NULL, outdim=NULL, method="box") {
            
            # Validate method parameter
            if (method != "box") {
              stop("Only 'box' method is currently supported")
            }
            
            # Validate input is 4D
            if (length(dim(x)) != 4) {
              stop("Input must be a 4D DenseNeuroVec object")
            }
            
            if (sum(!is.null(spacing), !is.null(factor), !is.null(outdim)) != 1) {
              stop("Exactly one of 'spacing', 'factor', or 'outdim' must be specified")
            }
            
            # Get current dimensions and spacing
            current_dims <- dim(x)
            current_space <- space(x)
            current_spacing <- current_space@spacing
            
            # Calculate new dimensions
            new_dims <- calculate_downsample_dims(current_dims, current_spacing,
                                                 spacing, factor, outdim)
            
            # Perform downsampling using C++ function
            downsampled_data <- downsample_4d_cpp(as.array(x), 
                                                 as.integer(new_dims),
                                                 as.integer(current_dims))
            
            # Calculate new spacing
            spatial_scale_factors <- current_dims[1:3] / new_dims[1:3]
            new_spacing <- current_spacing[1:3] * spatial_scale_factors
            
            # Add back the time spacing if 4D
            if (length(current_spacing) == 4) {
              new_spacing <- c(new_spacing, current_spacing[4])
            }
            
            # Create new NeuroSpace with updated dimensions and spacing
            new_space <- NeuroSpace(dim = new_dims,
                                  spacing = new_spacing,
                                  origin = current_space@origin,
                                  axes = current_space@axes,
                                  trans = current_space@trans)
            
            # Return new DenseNeuroVec
            DenseNeuroVec(downsampled_data, new_space, label=x@label)
          })

#' Downsample a NeuroVec
#'
#' @rdname downsample-methods
#' @export
setMethod(f="downsample", signature=signature("NeuroVec"),
          def=function(x, spacing=NULL, factor=NULL, outdim=NULL, method="box") {
            # For generic NeuroVec, try to convert to DenseNeuroVec first
            if (is(x, "SparseNeuroVec")) {
              stop("Downsampling of SparseNeuroVec not yet implemented. Convert to DenseNeuroVec first.")
            }
            callGeneric(x, spacing, factor, outdim, method)
          })