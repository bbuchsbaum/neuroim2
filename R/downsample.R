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

#' @keywords internal
#' @noRd
.downsample_bin_axis <- function(coord, old_n, new_n) {
  floor((coord - 1) * new_n / old_n) + 1L
}

#' @keywords internal
#' @noRd
.downsample_sparse_coords <- function(coords, old_dims, new_dims) {
  cbind(
    .downsample_bin_axis(coords[, 1], old_dims[1], new_dims[1]),
    .downsample_bin_axis(coords[, 2], old_dims[2], new_dims[2]),
    .downsample_bin_axis(coords[, 3], old_dims[3], new_dims[3])
  )
}

#' @keywords internal
#' @noRd
.downsample_linear_index_3d <- function(coords, dims) {
  1L +
    (coords[, 1] - 1L) +
    (coords[, 2] - 1L) * dims[1] +
    (coords[, 3] - 1L) * dims[1] * dims[2]
}

#' @keywords internal
#' @noRd
.downsample_sparse_neurovec <- function(x, new_dims) {
  current_dims <- dim(x)
  spatial_dims <- current_dims[1:3]
  nt <- current_dims[4]
  input_idx <- indices(x)
  input_coords <- arrayInd(input_idx, .dim = spatial_dims)

  if (nrow(input_coords) == 0L) {
    stop("Cannot downsample a SparseNeuroVec with an empty mask.")
  }

  output_coords <- .downsample_sparse_coords(input_coords, spatial_dims, new_dims[1:3])
  output_keys <- do.call(paste, c(as.data.frame(output_coords), sep = ":"))
  voxel_groups <- split(seq_len(nrow(output_coords)), output_keys)

  first_members <- vapply(voxel_groups, `[`, integer(1), 1L)
  kept_coords <- output_coords[first_members, , drop = FALSE]
  kept_linear_idx <- .downsample_linear_index_3d(kept_coords, new_dims[1:3])
  order_idx <- order(kept_linear_idx)
  voxel_groups <- voxel_groups[order_idx]
  kept_coords <- kept_coords[order_idx, , drop = FALSE]
  kept_linear_idx <- kept_linear_idx[order_idx]

  out_data <- vapply(
    voxel_groups,
    function(idx) {
      rowMeans(x@data[, idx, drop = FALSE])
    },
    numeric(nt)
  )

  if (!is.matrix(out_data)) {
    out_data <- matrix(out_data, nrow = nt, ncol = 1L)
  }
  colnames(out_data) <- NULL

  out_mask <- array(FALSE, dim = new_dims[1:3])
  out_mask[kept_linear_idx] <- TRUE

  current_space <- space(x)
  current_spacing <- current_space@spacing
  spatial_scale_factors <- current_dims[1:3] / new_dims[1:3]
  new_spacing <- current_spacing[1:3] * spatial_scale_factors
  new_trans <- rescale_affine(
    current_space@trans,
    shape = current_dims[1:3],
    zooms = new_spacing,
    new_shape = new_dims[1:3]
  )

  if (length(current_spacing) == 4) {
    new_spacing <- c(new_spacing, current_spacing[4])
  }

  new_space <- NeuroSpace(
    dim = new_dims,
    spacing = new_spacing,
    origin = new_trans[1:3, 4],
    axes = current_space@axes,
    trans = new_trans
  )

  SparseNeuroVec(out_data, new_space, out_mask, label = x@label)
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
            new_trans <- rescale_affine(
              current_space@trans,
              shape = current_dims[1:3],
              zooms = new_spacing,
              new_shape = new_dims[1:3]
            )
            
            # Add back the time spacing if 4D
            if (length(current_spacing) == 4) {
              new_spacing <- c(new_spacing, current_spacing[4])
            }
            
            # Create new NeuroSpace with updated dimensions and spacing
            new_space <- NeuroSpace(dim = new_dims,
                                  spacing = new_spacing,
                                  origin = new_trans[1:3, 4],
                                  axes = current_space@axes,
                                  trans = new_trans)
            
            # Return new DenseNeuroVec
            DenseNeuroVec(downsampled_data, new_space, label=x@label)
          })

#' Downsample a SparseNeuroVec
#'
#' @rdname downsample-methods
#' @export
setMethod(f = "downsample", signature = signature(x = "SparseNeuroVec"),
          def = function(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box") {
            if (method != "box") {
              stop("Only 'box' method is currently supported")
            }

            if (length(dim(x)) != 4) {
              stop("Input must be a 4D SparseNeuroVec object")
            }

            if (sum(!is.null(spacing), !is.null(factor), !is.null(outdim)) != 1) {
              stop("Exactly one of 'spacing', 'factor', or 'outdim' must be specified")
            }

            current_dims <- dim(x)
            current_spacing <- space(x)@spacing
            new_dims <- calculate_downsample_dims(
              current_dims,
              current_spacing,
              spacing,
              factor,
              outdim
            )

            .downsample_sparse_neurovec(x, new_dims)
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

#' Calculate new dimensions for 3D downsampling
#' @keywords internal
#' @noRd
calculate_downsample_dims_3d <- function(current_dims, current_spacing, 
                                        spacing = NULL, factor = NULL, outdim = NULL) {
  
  # current_dims should be 3D (x, y, z)
  spatial_dims <- current_dims
  
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
  
  # Return 3D dimensions
  new_spatial_dims
}

#' Downsample a DenseNeuroVol
#'
#' @rdname downsample-methods
#' @param x A DenseNeuroVol object to downsample
#' @param spacing Target voxel spacing (numeric vector of length 3)
#' @param factor Downsampling factor (single value or vector of length 3, between 0 and 1)
#' @param outdim Target output dimensions (numeric vector of length 3)
#' @param method Downsampling method (currently only "box" for box averaging)
#'
#' @examples
#' # Create a sample 3D volume
#' data <- array(rnorm(64*64*32), dim = c(64, 64, 32))
#' space <- NeuroSpace(dim = c(64, 64, 32), 
#'                     origin = c(0, 0, 0),
#'                     spacing = c(2, 2, 2))
#' vol <- DenseNeuroVol(data, space)
#'
#' # Downsample by factor
#' vol_down1 <- downsample(vol, factor = 0.5)
#' 
#' # Downsample to target spacing
#' vol_down2 <- downsample(vol, spacing = c(4, 4, 4))
#' 
#' # Downsample to target dimensions
#' vol_down3 <- downsample(vol, outdim = c(32, 32, 16))
#'
#' @export
setMethod(f="downsample", signature=signature("DenseNeuroVol"),
          def=function(x, spacing=NULL, factor=NULL, outdim=NULL, method="box") {
            
            # Validate method parameter
            if (method != "box") {
              stop("Only 'box' method is currently supported")
            }
            
            # Validate input is 3D
            if (length(dim(x)) != 3) {
              stop("Input must be a 3D DenseNeuroVol object")
            }
            
            if (sum(!is.null(spacing), !is.null(factor), !is.null(outdim)) != 1) {
              stop("Exactly one of 'spacing', 'factor', or 'outdim' must be specified")
            }
            
            # Get current dimensions and spacing
            current_dims <- dim(x)
            current_space <- space(x)
            current_spacing <- current_space@spacing[1:3]  # Take only spatial spacing
            
            # Calculate new dimensions
            new_dims <- calculate_downsample_dims_3d(current_dims, current_spacing,
                                                    spacing, factor, outdim)
            
            # Perform downsampling using C++ function
            downsampled_data <- downsample_3d_cpp(as.array(x), 
                                                 as.integer(new_dims),
                                                 as.integer(current_dims))
            
            # Calculate new spacing
            spatial_scale_factors <- current_dims / new_dims
            new_spacing <- current_spacing * spatial_scale_factors
            new_trans <- rescale_affine(
              current_space@trans,
              shape = current_dims,
              zooms = new_spacing,
              new_shape = new_dims
            )
            
            # Create new NeuroSpace with updated dimensions and spacing
            new_space <- NeuroSpace(dim = new_dims,
                                  spacing = new_spacing,
                                  origin = new_trans[1:3, 4],
                                  axes = current_space@axes,
                                  trans = new_trans)
            
            # Return new DenseNeuroVol
            DenseNeuroVol(downsampled_data, new_space, label=x@label)
          })

#' Downsample a NeuroVol
#'
#' @rdname downsample-methods
#' @export
setMethod(f="downsample", signature=signature("NeuroVol"),
          def=function(x, spacing=NULL, factor=NULL, outdim=NULL, method="box") {
            # For generic NeuroVol, try to convert to DenseNeuroVol first
            if (is(x, "SparseNeuroVol")) {
              stop("Downsampling of SparseNeuroVol not yet implemented. Convert to DenseNeuroVol first.")
            }
            callGeneric(x, spacing, factor, outdim, method)
          })
