#' Five-Dimensional Sparse Neuroimaging Data Container
#'
#' @title NeuroHyperVec Class for Multi-Feature Neuroimaging Data
#' @description
#' The \code{NeuroHyperVec} class provides an efficient container for five-dimensional
#' neuroimaging data where spatial dimensions are sparse. It is particularly suited for
#' analyses involving multiple features per trial/timepoint, such as basis functions,
#' spectral components, or multi-modal measurements.
#'
#' @slot mask A \code{\linkS4class{LogicalNeuroVol}} object defining the spatial mask.
#' @slot data A three-dimensional array with dimensions [features x trials x voxels] containing the data.
#' @slot space A \code{\linkS4class{NeuroSpace}} object defining the 5D space.
#' @slot lookup_map An integer vector for O(1) spatial index lookups.
#'
#' @details
#' The class organizes data in a 5D structure:
#' \itemize{
#'   \item Dimensions 1-3: Spatial coordinates (x, y, z)
#'   \item Dimension 4: Trials or timepoints
#'   \item Dimension 5: Features or measurements
#' }
#'
#' Data is stored internally as a three-dimensional array for efficiency:
#' \itemize{
#'   \item Dimensions 1: Features (dimension 5)
#'   \item Dimensions 2: Trials (dimension 4)
#'   \item Dimensions 3: Voxels (flattened spatial)
#' }
#'
#' Key features:
#' \itemize{
#'   \item Memory-efficient sparse storage of spatial dimensions
#'   \item Fast access to feature vectors and time series
#'   \item Flexible indexing across all dimensions
#'   \item Maintains spatial relationships and metadata
#' }
#'
#' @seealso
#' \code{\linkS4class{NeuroVec}}, \code{\linkS4class{LogicalNeuroVol}}, \code{\linkS4class{NeuroSpace}}
#'
#' @examples
#'
#' # Create a simple 5D dataset (10x10x10 spatial, 5 trials, 3 features)
#' dims <- c(10, 10, 10)
#' space <- NeuroSpace(c(dims, 5, 3))
#'
#' # Create a sparse mask (20% of voxels)
#' mask_data <- array(runif(prod(dims)) < 0.2, dims)
#' mask <- LogicalNeuroVol(mask_data, NeuroSpace(dims))
#'
#' # Generate random data for active voxels
#' n_voxels <- sum(mask_data)
#' data <- array(rnorm(3 * 5 * n_voxels), dim = c(3, 5, n_voxels))  # [features x trials x voxels]
#'
#' # Create NeuroHyperVec object
#' hvec <- NeuroHyperVec(data, space, mask)
#'
#' # Access operations
#' # Get data for specific voxel across all trials/features
#' series(hvec, 5, 5, 5)
#'
#' # Extract a 3D volume for specific trial and feature
#' hvec[,,,2,1]
#'
#'
#' @name NeuroHyperVec-class
#' @rdname NeuroHyperVec-class
#' @exportClass NeuroHyperVec
setClass(
  "NeuroHyperVec",
  slots = c(
    mask = "LogicalNeuroVol",
    data = "array",  # Dimensions: [features x trials x voxels]
    space = "NeuroSpace",
    lookup_map = "integer"
  ),
  contains = c( "ArrayLike5D"),
  validity = function(object) {
    # Validate that data is a 3D array
    if (!is.array(object@data) || length(dim(object@data)) != 3) {
      return("Data must be a 3D array with dimensions [features x trials x voxels]")
    }

    # Get expected dimensions
    num_voxels <- sum(object@mask@.Data)
    num_trials <- dim(object@space)[4]
    num_features <- dim(object@space)[5]

    # Validate array dimensions
    expected_dims <- c(num_features, num_trials, num_voxels)
    if (!identical(dim(object@data), expected_dims)) {
      return(sprintf(
        "Data array dimensions [%s] do not match expected [%d x %d x %d] (features x trials x voxels)",
        paste(dim(object@data), collapse=" x "),
        num_features, num_trials, num_voxels))
    }

    # Validate that mask is consistent with spatial dimensions
    mask_dims <- dim(object@mask)
    if (length(mask_dims) != 3) {
      return("Mask must be a 3D volume")
    }
    if (!all(mask_dims == dim(object@space)[1:3])) {
      return("Mask dimensions must match spatial dimensions of the space")
    }

    # Validate lookup_map
    if (length(object@lookup_map) != prod(dim(object@mask))) {
      return("lookup_map length must match the total number of voxels in the mask")
    }

    TRUE
  }
)

#' Constructor for NeuroHyperVec class
#'
#' @param data A matrix or three-dimensional array containing the data.
#' @param space A \code{\linkS4class{NeuroSpace}} object defining the spatial dimensions.
#' @param mask A mask volume (array, vector, or \code{\linkS4class{LogicalNeuroVol}}).
#' @return A new \code{\linkS4class{NeuroHyperVec}} object.
#'
#' @seealso \code{\linkS4class{NeuroSpace}}, \code{\linkS4class{LogicalNeuroVol}}
#'
#' @examples
#' # Create a 5D space (10x10x10 spatial, 2 trials, 2 features)
#' space <- NeuroSpace(c(10,10,10,2,2))
#'
#' # Create a mask for the spatial dimensions
#' space3d <- NeuroSpace(c(10,10,10))
#' mask_data <- array(TRUE, dim=c(10,10,10))  # All voxels active
#' mask <- LogicalNeuroVol(mask_data, space3d)
#'
#' # Create data in the format expected by NeuroHyperVec:
#' # 3D array with dimensions [features x trials x voxels]
#' n_features <- 2
#' n_trials <- 2
#' n_voxels <- sum(mask_data)  # 1000 voxels
#' data_array <- array(rnorm(n_features * n_trials * n_voxels),
#'                    dim = c(n_features, n_trials, n_voxels))
#'
#' # Create the NeuroHyperVec object
#' hvec <- NeuroHyperVec(data_array, space, mask)
#'
#' @export
NeuroHyperVec <- function(data, space, mask) {
  # Validate input types
  if (!is.matrix(data) && !is.array(data)) {
    stop("'data' must be a matrix or a three-dimensional array.")
  }
  if (!inherits(space, "NeuroSpace")) {
    stop("'space' must be a NeuroSpace object")
  }

  # Get dimensions from space
  dims <- dim(space)
  num_features <- dims[5]
  num_trials <- dims[4]

  # Convert mask to LogicalNeuroVol if needed
  if (is.array(mask)) {
    if (length(dim(mask)) != 3) {
      stop("'mask' must be a 3D logical array")
    }
    mask_data <- as.vector(mask)
    mask <- LogicalNeuroVol(mask, NeuroSpace(dim(mask)))
  } else if (is.vector(mask) && is.logical(mask)) {
    mask_space <- NeuroSpace(dim = dims[1:3])
    mask <- LogicalNeuroVol(mask, mask_space)
  } else if (!inherits(mask, "LogicalNeuroVol")) {
    stop("'mask' must be a 3D logical array, a 1D logical vector, or an instance of 'LogicalNeuroVol'.")
  }

  # Validate mask dimensions
  if (!all(dim(mask)[1:3] == dims[1:3])) {
    stop("Mask dimensions must match spatial dimensions of the space")
  }

  # Get number of voxels
  spatial_dims <- dims[1:3]
  num_spatial_voxels <- prod(spatial_dims)
  num_mask_voxels <- sum(mask@.Data)

  # Validate data dimensions
  if (is.array(data) && !is.matrix(data)) {
    data_dims <- dim(data)
    expected_dims <- c(num_features, num_trials, num_mask_voxels)
    if (!all(dim(data) == expected_dims)) {
      stop(sprintf("Data array dimensions [%s] do not match expected [%s]",
                  paste(data_dims, collapse = " x "),
                  paste(expected_dims, collapse = " x ")))
    }
  }

  # Handle matrix input
  if (is.matrix(data)) {
    if (ncol(data) == num_mask_voxels) {
      # Case: features x voxels or (features * trials) x voxels
      if (nrow(data) == num_features) {
        # Case: features x voxels
        data_array <- array(0, dim=c(num_features, num_trials, num_mask_voxels))
        for (t in 1:num_trials) {
          data_array[,t,] <- data
        }
      } else if (nrow(data) == num_features * num_trials) {
        # Case: (features * trials) x voxels
        data_array <- array(0, dim=c(num_features, num_trials, num_mask_voxels))
        for (t in 1:num_trials) {
          idx_start <- (t-1) * num_features + 1
          idx_end <- t * num_features
          data_array[,t,] <- data[idx_start:idx_end,]
        }
      } else {
        stop("For matrix input with ncol == num_mask_voxels, nrow must be either num_features or (num_features * num_trials)")
      }
    } else if (ncol(data) == num_mask_voxels * num_trials) {
      # Case: features x (trials * voxels)
      if (nrow(data) != num_features) {
        stop("For matrix input with ncol == num_mask_voxels * num_trials, nrow must be num_features")
      }
      data_array <- array(0, dim=c(num_features, num_trials, num_mask_voxels))
      for (t in 1:num_trials) {
        idx_start <- (t-1) * num_mask_voxels + 1
        idx_end <- t * num_mask_voxels
        data_array[,t,] <- data[,idx_start:idx_end]
      }
    } else {
      stop("Invalid matrix dimensions. Must be either:
           - features x voxels
           - features x (trials * voxels)
           - (features * trials) x voxels")
    }
  } else {
    # Handle array input
    if (length(dim(data)) != 3) {
      stop("'data' array must have exactly three dimensions [features x trials x voxels]")
    }
    data_array <- data
    if (dim(data)[1] != num_features) {
      stop("First dimension of 'data' must match number of features")
    }
    if (dim(data)[2] != num_trials) {
      stop("Second dimension of 'data' must match number of trials")
    }
    if (dim(data)[3] != num_mask_voxels) {
      stop("Third dimension of 'data' must match the number of non-zero mask elements")
    }
  }

  # Create lookup map
  num_spatial_voxels <- prod(dim(space)[1:3])
  lookup_map <- integer(num_spatial_voxels)
  lookup_map[which(mask@.Data)] <- seq_len(num_mask_voxels)

  # Create and return new object
  new("NeuroHyperVec",
      data = data_array,
      space = space,
      mask = mask,
      lookup_map = lookup_map)
}

#' Series method for NeuroHyperVec
#'
#' @param x The NeuroHyperVec object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... Additional arguments (not used)
#'
#' @details when x is a NeuroHyperVec object, the series method returns a 2D array with dimensions [features x trials]
#' @return A 2D array with dimensions [features x trials]
setMethod("series", signature(x = "NeuroHyperVec"),
  function(x, i, j, k, ...) {
    # Validate indices
    if (missing(i) || missing(j) || missing(k)) {
      stop("Must provide indices 'i', 'j', and 'k'")
    }
    spatial_dims <- dim(x@space)[1:3]
    if (i < 1 || i > spatial_dims[1] || j < 1 || j > spatial_dims[2] || k < 1 || k > spatial_dims[3]) {
      stop("Indices out of bounds")
    }
    # Convert 3D indices to linear indices
    linear_idx <- (k - 1) * (spatial_dims[1] * spatial_dims[2]) + (j - 1) * spatial_dims[1] + i

    # Get lookup index
    lookup_index <- x@lookup_map[linear_idx]

    # Check if voxel is in the mask
    if (lookup_index == 0) {
      return(array(0, dim = c(dim(x@data)[1], dim(x@data)[2])))
    }

    # Extract data
    data_slice <- x@data[,,lookup_index]

    data_slice
  }
)

#' Linear access method for NeuroHyperVec
#'
#' @param x The NeuroHyperVec object
#' @param i The linear indices
#' @param ... Additional arguments (not used)
#'
#' @rdname linear_access-methods
setMethod("linear_access", signature(x = "NeuroHyperVec"),
  function(x, i, ...) {
    # Get dimensions
    dims <- dim(x@space)
    num_spatial_voxels <- prod(dims[1:3])
    num_trials <- dims[4]
    num_features <- dims[5]
    total_elements <- num_spatial_voxels * num_trials * num_features

    # Validate indices
    if (any(i < 1) || any(i > total_elements)) {
      stop("indices must be within range of data dimensions")
    }

    tmp <- i - 1
    spatial_idx <- (tmp %% num_spatial_voxels) + 1
    tmp <- tmp %/% num_spatial_voxels
    trial_idx <- (tmp %% num_trials) + 1
    feature_idx <- (tmp %/% num_trials) + 1

    # Get lookup indices
    lookup_indices <- x@lookup_map[spatial_idx]

    # Create output
    out <- numeric(length(i))

    # Fill non-zero entries
    non_zero <- lookup_indices > 0
    if (any(non_zero)) {
      out[non_zero] <- x@data[
        cbind(
          feature_idx[non_zero],
          trial_idx[non_zero],
          lookup_indices[non_zero]
        )
      ]
    }

    out
  }
)

#' [ method for NeuroHyperVec
#'
#' @param x The NeuroHyperVec object
#' @param i,j,k,l,m Indices for each dimension
#' @param drop Whether to drop dimensions of length 1
#' @param ... Additional arguments (not used)
#'
#' @rdname NeuroHyperVec-class
#' @aliases [.NeuroHyperVec
setMethod("[", signature(x = "NeuroHyperVec"),
  function(x, i, j, k, l, m, ..., drop = TRUE) {
    dims <- dim(x@space)

    # Handle missing indices
    if (missing(i)) i <- seq_len(dims[1])
    if (missing(j)) j <- seq_len(dims[2])
    if (missing(k)) k <- seq_len(dims[3])
    if (missing(l)) l <- seq_len(dims[4])
    if (missing(m)) m <- seq_len(dims[5])

    # Validate indices
    validate_indices <- function(idx, max_val, name) {
      if (any(idx < 1) || any(idx > max_val)) {
        stop(sprintf("Index '%s' out of bounds", name))
      }
    }

    validate_indices(i, dims[1], "i")
    validate_indices(j, dims[2], "j")
    validate_indices(k, dims[3], "k")
    validate_indices(l, dims[4], "l")
    validate_indices(m, dims[5], "m")

    # Calculate spatial indices
    spatial_dims <- dims[1:3]
    spatial_idx_grid <- expand.grid(i = i, j = j, k = k)
    spatial_linear_idx <- with(spatial_idx_grid,
                               (k - 1) * (spatial_dims[1] * spatial_dims[2]) +
                               (j - 1) * spatial_dims[1] + i)

    # Get lookup indices
    lookup_indices <- x@lookup_map[spatial_linear_idx]

    # Initialize output array with correct dimensions
    out_dims <- if (drop) {
      c(length(i), length(j), length(k))
    } else {
      c(length(i), length(j), length(k), length(l), length(m))
    }
    out_array <- array(0, dim = out_dims)

    # Fill non-zero entries
    non_zero <- lookup_indices > 0
    if (any(non_zero)) {
      # Get positions where data should be filled
      non_zero_positions <- which(non_zero)

      # Map to multi-dimensional indices
      array_indices <- arrayInd(non_zero_positions, .dim = c(length(i), length(j), length(k)))

      # Retrieve data from x@data
      voxel_indices <- lookup_indices[non_zero]

      # For each non-zero position
      for (idx in seq_along(voxel_indices)) {
        ii <- array_indices[idx, 1]
        jj <- array_indices[idx, 2]
        kk <- array_indices[idx, 3]

        # Get the data for this voxel
        voxel_data <- x@data[m, l, voxel_indices[idx], drop = FALSE]

        if (drop) {
          out_array[ii, jj, kk] <- voxel_data[1, 1, 1]
        } else {
          out_array[ii, jj, kk, , ] <- voxel_data[, , 1]
        }
      }
    }

    out_array
  }
)


#' @rdname show-methods
#' @export
setMethod("show", signature(object="NeuroHyperVec"),
          def=function(object) {
            cat("\n", crayon::bold(crayon::blue("NeuroHyperVec Object")), "\n")
            cat(crayon::silver("======================================\n"))

            # Dimensions section
            dims <- dim(object@space)
            spatial_dims <- paste(dims[1:3], collapse=" x ")
            cat("\n", crayon::yellow("Dimensions:"), "\n")
            cat(" ", crayon::silver("."), " Spatial: ",
                crayon::green(spatial_dims),
                crayon::silver(" (xyz)"), "\n")
            cat(" ", crayon::silver("."), " Trials:  ",
                crayon::green(sprintf("%-6d", dims[4])),
                crayon::silver(" (4th dimension)"), "\n")
            cat(" ", crayon::silver("."), " Features:",
                crayon::green(sprintf("%-6d", dims[5])),
                crayon::silver(" (5th dimension)"), "\n")

            # Sparsity information
            n_total <- prod(dims[1:3])
            n_active <- sum(object@mask@.Data)
            sparsity <- round(100 * n_active / n_total, 2)
            cat("\n", crayon::yellow("Sparsity:"), "\n")
            cat(" ", crayon::silver("."), " Active Voxels: ",
                crayon::green(sprintf("%d / %d", n_active, n_total)), "\n")
            cat(" ", crayon::silver("."), " Coverage:      ",
                crayon::green(sprintf("%.2f%%", sparsity)), "\n")
            cat(" ", crayon::silver("."), " Compression:   ",
                crayon::green(sprintf("%.2fx", n_total/n_active)), "\n")

            # Memory usage
            data_size <- object.size(object@data)
            total_size <- object.size(object)
            cat("\n", crayon::yellow("Memory Usage:"), "\n")
            cat(" ", crayon::silver("."), " Data:    ",
                crayon::green(format(data_size, units="auto")), "\n")
            cat(" ", crayon::silver("."), " Total:   ",
                crayon::green(format(total_size, units="auto")), "\n")
            cat(" ", crayon::silver("."), " Overhead:",
                crayon::green(format(total_size - data_size, units="auto")), "\n")

            # Space information
            sp <- object@space
            spacing <- spacing(sp)[1:3]
            origin <- origin(sp)[1:3]
            cat("\n", crayon::yellow("Spatial Information:"), "\n")
            cat(" ", crayon::silver("."), " Spacing: ",
                crayon::green(sprintf("%.2f x %.2f x %.2f", spacing[1], spacing[2], spacing[3])),
                crayon::silver(" mm"), "\n")
            cat(" ", crayon::silver("."), " Origin:  ",
                crayon::green(sprintf("%.1f, %.1f, %.1f", origin[1], origin[2], origin[3])), "\n")

            # Footer with usage hints
            cat(crayon::silver("\n======================================\n"))
            cat("\n", crayon::bold("Access Methods:"), "\n")
            cat(" ", crayon::silver("."), " Extract Series:  ",
                crayon::blue("series(object, i, j, k)"), "\n")
            cat(" ", crayon::silver("."), " Get Subset:     ",
                crayon::blue("object[i, j, k, trial, feature]"), "\n")
            cat(" ", crayon::silver("."), " Get Volume:     ",
                crayon::blue("object[, , , 1, 1]"),
                crayon::silver(" # first trial, first feature"), "\n")
            cat(" ", crayon::silver("."), " Get Timeseries: ",
                crayon::blue("series(object, 10, 20, 30)"),
                crayon::silver(" # at xyz=(10,20,30)"), "\n\n")
          })

#' @rdname mask-methods
#' @export
setMethod("mask", "NeuroHyperVec",
          function(x) {
            x@mask
          })
