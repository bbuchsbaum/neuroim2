#' @include all_class.R
NULL
#' @include all_generic.R
NULL



#' SparseNeuroVecSource
#'
#' Constructs a SparseNeuroVecSource object
#'
#' @param meta_info an object of class \code{\linkS4class{MetaInfo}}
#' @param indices an optional vector of 1D indices the subset of volumes to load
#' @param mask a logical 3D \code{array},  a logical 1D \code{vector} or a \code{LogicalNeuroVol}
#' @return A \code{SparseNeuroVecSource} object
#' @keywords internal
#' @noRd
SparseNeuroVecSource <- function(meta_info, indices=NULL, mask) {

  # Input validation
  assert_that(length(dim(meta_info)) >= 3,
              msg = "meta_info must have at least 3 dimensions")

  indices <- if(is.null(indices)) seq(1, dim(meta_info)[4]) else indices
  assert_that(all(indices >= 1 & indices <= dim(meta_info)[4]),
              msg = "indices must be within valid range")

	D <- dim(meta_info)[1:3]


	if (is.vector(mask) && length(mask) < prod(D)) {
    ### this is a vector of indices
		m <- array(FALSE, D)
		m[mask] <- TRUE
		mask <- m
	} else if (identical(dim(mask), as.integer(D))) {
		mask <- as.array(mask)
	} else if (is.vector(mask) && length(mask) == prod(D)) {
		mask <- array(mask, D)
	} else {
		stop(sprintf("Invalid mask dimensions: %s", paste(dim(mask), collapse=" x ")))
	}

  if (!inherits(mask, "LogicalNeuroVol")) {
    mspace <- NeuroSpace(dim(mask),  meta_info@spacing, meta_info@origin, meta_info@spatial_axes)
    mask <- LogicalNeuroVol(mask, mspace)
  }

	assert_that(all(dim(mask) == D),
              msg = "mask dimensions must match data dimensions")

	new("SparseNeuroVecSource", meta_info=meta_info, indices=indices, mask=mask)
}

#' @noRd
prep_sparsenvec <- function(data, space, mask) {
  if (!inherits(mask, "LogicalNeuroVol")) {
    mspace <- NeuroSpace(dim(space)[1:3],
                         spacing(space),
                         origin(space),
                         axes(space),
                         trans(space))
    mask <- LogicalNeuroVol(as.logical(mask), mspace)
  }

  cardinality <- sum(mask)

  assert_that(inherits(mask, "LogicalNeuroVol"),
              msg = "mask must be a LogicalNeuroVol object")

  if (is.matrix(data)) {
    Nind <- sum(mask == TRUE)
    if (nrow(data) == Nind) {
      data <- t(data)
      assert_that(ncol(data) == cardinality, msg = "data matrix must match cardinality of `mask`")
    } else if (ncol(data) == Nind) {
      assert_that(ncol(data) == cardinality, msg = "data matrix must match cardinality of `mask`")
    } else {
      stop(sprintf("Matrix dimensions %s do not match mask cardinality %d",
                   paste(dim(data), collapse=" x "), Nind))
    }
    D4 <- dim(data)[1]
  } else if (length(dim(data)) == 4) {
    dims <- dim(data)
    data_mat <- matrix(0, nrow = dims[4], ncol = prod(dims[1:3]))
    for (t in 1:dims[4]) {
      vol_t <- data[,,,t]
      data_mat[t,] <- as.vector(vol_t)
    }
    data <- data_mat[, mask@.Data, drop=FALSE]  # Only keep masked voxels
    D4 <- dims[4]
  } else {
    stop("Data must be either a matrix or a 4D array.")
  }

  if (ndim(space) == 3) {
    space <- add_dim(space, D4)
  }

  stopifnot(ndim(space) == 4)

  list(mask = mask, data = data, space = space)
}

#' Construct a SparseNeuroVec Object
#'
#' Constructs a SparseNeuroVec object for efficient representation and manipulation
#' of sparse neuroimaging data with many zero or missing values.
#'
#' @param data A matrix or a 4-D array containing the neuroimaging data. The dimensions of the data should be consistent with the dimensions of the provided NeuroSpace object and mask.
#' @param space A \link{NeuroSpace} object representing the dimensions and voxel spacing of the neuroimaging data.
#' @param mask A 3D array, 1D vector of type logical, or an instance of type \link{LogicalNeuroVol}, which specifies the locations of the non-zero values in the data.
#' @param label Optional character string providing a label for the vector
#' @return A SparseNeuroVec object, containing the sparse neuroimaging data, mask, and associated NeuroSpace information.
#' @export
#'
#' @examples
#' bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
#' mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
#' mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
#' svec <- SparseNeuroVec(mat, bspace, mask)
#' length(indices(svec)) == sum(mask)
#' @rdname SparseNeuroVec-class
SparseNeuroVec <- function(data, space, mask, label = "") {
	stopifnot(inherits(space, "NeuroSpace"))
  
  # Ensure space has 4 dimensions
  if (ndim(space) != 4) {
    stop("The 'space' argument must have exactly 4 dimensions")
  }
  
  p <- prep_sparsenvec(data, space, mask)

	new("SparseNeuroVec", space=p$space, mask=p$mask,
	    map=IndexLookupVol(space(p$mask), as.integer(which(p$mask))), data=p$data, label=label)
}

#' @rdname load_data-methods
#' @export
setMethod(f="load_data", signature=c("SparseNeuroVecSource"),
		def=function(x) {

			meta <- x@meta_info
			nels <- prod(dim(meta)[1:3])

			ind <- x@indices
			M <- x@mask > 0
			reader <- data_reader(meta, offset=0)
			dat4D <- read_elements(reader, prod(dim(meta)[1:4]))
			close(reader)

			datlist <- lapply(1:length(ind), function(i) {
				offset <- (nels * (ind[i]-1))
				dat4D[(offset+1):(offset + nels)][M]
			})

			#close(reader)
			arr <- do.call(rbind, datlist)

			if (.hasSlot(meta, "slope")) {
			  if (meta@slope != 0) {
			    arr <- arr*meta@slope
			  }
			}

			bspace <- NeuroSpace(c(dim(meta)[1:3], length(ind)), meta@spacing,
			                     meta@origin, meta@spatial_axes, trans=trans(meta))

			SparseNeuroVec(arr, bspace, x@mask)

		})


#' @rdname indices-methods
#' @keywords internal
setMethod(f="indices", signature=signature(x="AbstractSparseNeuroVec"),
          def=function(x) {
            indices(x@map)
          })


#' @export
#' @rdname coords-methods
setMethod(f="coords", signature=signature(x="AbstractSparseNeuroVec"),
          def=function(x,i) {
            if (missing(i)) {
              return(coords(x@map, indices(x@map)))
            }
            coords(x@map, i)
          })



#' @rdname series-methods
#' @export
setMethod("series", signature(x="AbstractSparseNeuroVec", i="ROICoords"),
          def=function(x,i) {
            callGeneric(x, coords(i))
          })


#' @export
#' @rdname series-methods
setMethod(f="series", signature=signature(x="AbstractSparseNeuroVec", i="matrix"),
         def=function(x,i) {
           idx <- grid_to_index(x@mask, i)
           callGeneric(x,idx)
         })


#' @export
#' @rdname series-methods
setMethod("series", signature(x="AbstractSparseNeuroVec", i="numeric"),
          def=function(x,i, j, k) {
            if (missing(j) && missing(k)) {
              callGeneric(x, as.integer(i))
            } else {
              callGeneric(x, as.integer(i), as.integer(j), as.integer(k))
            }
          })

#' @rdname series-methods
#' @export
setMethod(
  "series",
  signature(x = "AbstractSparseNeuroVec", i = "integer"),
  function(x, i, j, k, drop = TRUE) {
    # Case 1: user provided only i (voxel linear indices)
    if (missing(j) && missing(k)) {
      # Map linear indices -> actual row in sparse matrix or 0 if none
      mapped_idx <- lookup(x, i)  # vector of the same length as i
      # Prepare output: #rows = time, #cols = length(i)
      out <- matrix(0, nrow = dim(x)[4], ncol = length(i))

      # Identify which of those voxel indices are actually non-zero
      nz <- which(mapped_idx > 0)
      if (length(nz) > 0) {
        # Access the non-zero columns from x@data
        # Because x@data is (time x voxels)
        # We want to fill the columns out[, nz] from x@data[, mapped_idx[nz]]
        out[, nz] <- matricized_access(x, mapped_idx[nz])
      }

      # If user says drop=TRUE and asked for a single voxel, drop down to vector
      if (drop && length(i) == 1) {
        out <- drop(out)  # becomes a vector
      }
      return(out)

    } else {
      # Case 2: user gave i,j,k (3D coordinates).
      # same approach: convert (i,j,k) to linear indices, then do the same logic
      if (length(i) == 1 && length(j) == 1 && length(k) == 1) {
        # single voxel coordinate => we fill a single time-series
        # 1) convert (i,j,k) -> linear
        idx <- .gridToIndex3D(dim(x)[1:3], matrix(c(i, j, k), nrow=1))
        mapped_idx <- lookup(x, idx)
        out <- rep(0, dim(x)[4])  # time vector
        if (mapped_idx > 0) {
          # fill from x@data
          out <- x@data[, mapped_idx]
        }
        if (drop) {
          # already a 1D vector, so nothing special
          return(out)
        } else {
          # return a 2D matrix [time x 1]
          return(matrix(out, nrow = length(out), ncol = 1))
        }
      } else {
        # multiple 3D coords. Return a matrix [time x #coords]
        # expand (i, j, k) to a set of voxel indices
        coords_mat <- cbind(i, j, k)
        lin_idx <- .gridToIndex3D(dim(x)[1:3], coords_mat)
        mapped_idx <- lookup(x, lin_idx)
        out <- matrix(0, nrow = dim(x)[4], ncol = nrow(coords_mat))
        nz <- which(mapped_idx > 0)
        if (length(nz) > 0) {
          out[, nz] <- matricized_access(x, mapped_idx[nz])
        }
        # If user requested drop=TRUE and exactly one voxel, drop dimension
        if (drop && nrow(coords_mat) == 1) {
          out <- drop(out)  # becomes a vector of length = time
        }
        return(out)
      }
    }
  }
)


#' @param nonzero only include nonzero vectors in output list
#' @export
#' @rdname vectors-methods
#' @examples
#'
#' file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
#' vec <- read_vec(file_name)
#' v <- vectors(vec)
#' mean(v[[1]])
setMethod(f="vectors", signature=signature(x="SparseNeuroVec", subset="missing"),
          def = function(x, nonzero=FALSE) {
            if (nonzero) {
              force(x)
              ind <- indices(x)
              f <- function(i) series(x, ind[i])
              #lis <- lapply(seq_along(ind), function(i) f)
              deflist::deflist(f, length(ind))
            } else {
              ind <- 1:prod(dim(x)[1:3])
              vox <- index_to_grid(x, ind)
              f <- function(i) series(x, vox[i,1], vox[i,2], vox[i,3])
              #lis <- map(ind, function(i) f)
              deflist::deflist(f, length(ind))
            }

          })





#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="AbstractSparseNeuroVec", y="missing"),
          def=function(x,y,...) {
            x
          })


#' @export
#' @rdname concat-methods
setMethod(f="concat", signature=signature(x="SparseNeuroVec", y="SparseNeuroVec"),
          def=function(x,y,...) {
            if (!all(indices(x) == indices(y))) {
              stop("cannot concatenate arguments with different index maps")
            }

            if (!all(dim(x)[1:3] == dim(y)[1:3])) {
              stop("cannot concatenate arguments with different spatial dimensions")
            }

            ndat <- rbind(x@data, y@data)
            d1 <- dim(x)
            d2 <- dim(y)

            rest <- list(...)


            if (length(rest) >= 1) {
              mat <- do.call(rbind, map(rest, ~ .@data))

              ndim <- c(d1[1:3], d1[4] + d2[4] + nrow(mat))
              ndat <- rbind(ndat, mat)
              nspace <- NeuroSpace(ndim, spacing(x@space),  origin(x@space), axes(x@space), trans(x@space))
              SparseNeuroVec(ndat, nspace, mask=x@mask)
            } else {
              ndim <- c(d1[1:3], d1[4] + d2[4])
              nspace <- NeuroSpace(ndim, spacing(x@space),  origin(x@space), axes(x@space), trans(x@space))
              SparseNeuroVec(ndat, nspace, mask=x@mask)
            }

          })



#' @rdname lookup-methods
#' @keywords internal
setMethod(f="lookup", signature=signature(x="AbstractSparseNeuroVec", i="numeric"),
         def=function(x,i) {
            lookup(x@map, i)
          })

#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "matrix"),
          def=function (x, i) {
            x@data[i]
          })

#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "integer"),
          def=function (x, i) {
            x@data[,i]
          })

#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "numeric"),
          def=function (x, i) {
            x@data[,i]
          })


#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "matrix"),
          def=function (x, i) {
            x@data[i]
          })

#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "integer"),
          def=function (x, i) {
            x@data[,i]
          })

#' @export
#' @rdname matricized_access-methods
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "numeric"),
          def=function (x, i) {
            x@data[,i]
          })



#' @export
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x = "AbstractSparseNeuroVec", i = "numeric"),
  def = function(x, i) {
    # -------------------------------
    # Input Validation
    # -------------------------------
    if (!is.numeric(i)) {
      stop("'i' must be a numeric vector.")
    }

    if (any(is.na(i))) {
      stop("'i' contains NA values, which are not allowed.")
    }

    if (any(i <= 0)) {
      stop("All indices in 'i' must be positive integers.")
    }

    if (any(i != floor(i))) {
      stop("All indices in 'i' must be integers.")
    }

    # -------------------------------
    # Dimension Retrieval and Validation
    # -------------------------------
    dims <- dim(x)
    if (is.null(dims) || length(dims) < 4) {
      stop("The object 'x' must have at least 4 dimensions.")
    }

    spatial_nels <- prod(dims[1:3])  # Number of elements in the first three dimensions
    num_timepoints <- dims[4]        # Fourth dimension (e.g., time)

    # Total number of elements in 'x'
    total_elements <- spatial_nels * num_timepoints
    if (any(i > total_elements)) {
      stop(sprintf("Indices in 'i' exceed the total number of elements (%d).", total_elements))
    }

    # -------------------------------
    # Mapping Linear Indices to 4D Coordinates
    # -------------------------------
    # Calculate timepoints and spatial_offsets using integer division and modulo
    timepoints <- ((i - 1) %/% spatial_nels) + 1
    spatial_offsets <- ((i - 1) %% spatial_nels) + 1

    # -------------------------------
    # Sparse Lookup
    # -------------------------------
    # Perform lookup to get mapping indices; assumes 'lookup' returns 0 for zeros
    lookup_values <- lookup(x, spatial_offsets)

    # Identify non-zero lookups
    non_zero_indices <- which(lookup_values > 0)

    # Early exit if all lookups are zero
    if (length(non_zero_indices) == 0) {
      return(rep(0, length(i)))  # All requested values are zero
    }

    # -------------------------------
    # Prepare Indices for Data Retrieval
    # -------------------------------
    # Extract corresponding timepoints and spatial indices for non-zero lookups
    data_indices <- lookup_values[non_zero_indices]  # Indices in the sparse data matrix

    # Create a two-column matrix for 'matricized_access'
    idx_matrix <- cbind(data_indices, timepoints[non_zero_indices])

    # -------------------------------
    # Retrieve Non-Zero Values
    # -------------------------------
    # Retrieve the non-zero values from the data matrix
    non_zero_values <- matricized_access(x, idx_matrix)

    # -------------------------------
    # Assemble Output Vector
    # -------------------------------
    # Initialize the output vector with zeros
    output_values <- numeric(length(i))

    # Assign the retrieved non-zero values to their respective positions
    output_values[non_zero_indices] <- non_zero_values

    return(output_values)
  }
)


#setMethod(
#  f = "[",
#  signature = signature(x = "AbstractSparseNeuroVec", i = "numeric", j = "numeric"),
#  def = function(x, i, j, k, m, ..., drop = TRUE) {
#    # -------------------------------
#    # Input Validation
#    # -------------------------------
#    dims <- dim(x)
#    if (missing(k)) {
#      k <- seq_len(dims[3])
#    }
#    if (missing(m)) {
#      m <- seq_len(dims[4])
#    }
#
#    # Generate all combinations of spatial indices
#    grid <- expand.grid(i = i, j = j, k = k)
#    ind <- .gridToIndex3D(dims[1:3], as.matrix(grid))
#
#    # Perform lookup to get data indices in the sparse representation
#    mapped <- lookup(x, ind)
#
#    # Identify non-zero mappings
#    non_zero_positions <- which(mapped > 0)
#    if (length(non_zero_positions) == 0) {
#      result <- array(0, dim = c(length(i), length(j), length(k), length(m)))
#      if (drop) {
#        return(drop(result))
#      } else {
#        return(result)
#      }
#    }
#
#    # Extract mapped indices and corresponding grid positions for non-zero entries
#    data_indices <- mapped[non_zero_positions]  # Indices in the sparse data matrix
#    grid_indices <- grid[non_zero_positions, ]  # Corresponding spatial indices
#
#    # Retrieve data for specified voxels and timepoints
#    # Data is stored as [time, voxels], where each column is a voxel's time series
#    data_values <- x@data[m, data_indices, drop = FALSE]  # dimensions: n_m x n_voxels
#
#    # Initialize output array
#    result <- array(0, dim = c(length(i), length(j), length(k), length(m)))
#
#    # Calculate output positions
#    i_pos <- match(grid_indices$i, i)
#    j_pos <- match(grid_indices$j, j)
#    k_pos <- match(grid_indices$k, k)
#
#    # Fill the output array
#    for (idx in seq_along(non_zero_positions)) {
#        result[i_pos[idx], j_pos[idx], k_pos[idx], ] <- data_values[, idx]
#    }
#
#    if (drop) {
#        return(drop(result))
#    } else {
#        return(result)
#    }
#  }
#)


#' Extractor Method for AbstractSparseNeuroVec
#'
#' @description
#' Extracts a subset of data from a sparse four-dimensional brain image based on provided indices.
#'
#' @param x An object of class \code{AbstractSparseNeuroVec}
#' @param i Numeric vector specifying the indices for the first dimension
#' @param j Numeric vector specifying the indices for the second dimension
#' @param k Numeric vector specifying the indices for the third dimension (optional)
#' @param m Numeric vector specifying the indices for the fourth dimension (optional)
#' @param ... Additional arguments passed to methods
#' @param drop Logical indicating whether to drop dimensions of length one (default: TRUE)
#'
#' @return An array containing the extracted subset
#' 
#' @export
setMethod(f="[", signature=signature(x = "AbstractSparseNeuroVec", i = "numeric", j = "numeric"),
          def = function (x, i, j, k, m, ..., drop = TRUE) {
            if (missing(k))
              k = 1:(dim(x)[3])
            if (missing(m)) {
              m <- 1:(dim(x)[4])
            }

            vmat <- as.matrix(expand.grid(i,j,k))
            ind <- .gridToIndex3D(dim(x)[1:3], vmat[,1:3,drop = FALSE])

            mapped <- lookup(x, ind)
            keep <- mapped > 0
            dimout <- c(length(i),length(j),length(k),length(m))

            if (sum(keep) == 0) {
              if (drop) {
                return(drop(array(0, dimout)))
              } else {
                return(array(0, dimout))
              }
            }


            egrid <- expand.grid(mapped[keep], m)
            indmat <- cbind(egrid[,2], egrid[,1])

            oval <- numeric(prod(dimout))

            ## TODO assumes x has @data member ...
            ##oval[rep(keep, length(m))] <- x@data[indmat]
            oval[rep(keep, length(m))] <- matricized_access(x, indmat)

            dim(oval) <- c(length(i),length(j),length(k),length(m))

            if (drop) {
              drop(oval)
            } else {
              oval
            }
})


#' Extract a sub-vector
#' @name sub_vector
#' @rdname sub_vector-methods
#' @aliases sub_vector,SparseNeuroVec,numeric-method
#' @export
setMethod(f="sub_vector", signature=signature(x="SparseNeuroVec", i="numeric"),
          def=function(x, i) {
            assertthat::assert_that(max(i) <= dim(x)[4])
            
            # Get the subset of data for the requested timepoints
            res <- x@data[i,, drop=FALSE]
            
            # Create new space with updated dimensions
            xs <- space(x)
            newdim <- c(dim(x)[1:3], length(i))
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), 
                               axes(xs), trans(xs))
            
            # Create new SparseNeuroVec with subset of data
            new("SparseNeuroVec", space=bspace, mask=x@mask,
                map=x@map, data=res, label=x@label)
          })

#' [[
#'
#' @rdname SparseNeuroVec-methods
#' @param x the object
#' @param i the volume index
#' @export
setMethod(f="[[", signature=signature(x="SparseNeuroVec", i="numeric"),
          def = function(x, i) {
            stopifnot(length(i) == 1)
            xs <- space(x)
            dat <- x@data[i,]
            newdim <- dim(xs)[1:3]
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            SparseNeuroVol(dat, bspace, indices=indices(x))
          })



#' @name as
#' @export
setAs(from="SparseNeuroVec", to="matrix",
		  function(from) {
		    ind <- indices(from)
		    out <- matrix(0, dim(from)[4], prod(dim(from)[1:3]))
		    out[, ind] <- from@data
		    t(out)
			  #from@data
		  })


#' @export
setAs(from="SparseNeuroVec", to="DenseNeuroVec",
      function(from) {
        mat <- as(from, "matrix")
        DenseNeuroVec(mat, space(from))
      })


#' Convert to Matrix
#'
#' @aliases as.matrix,SparseNeuroVec-method
#'          as.matrix,NeuroVec-method
#' @param x The object to convert to a matrix
#' @param ... Additional arguments
#' @return A matrix representation of the object
#' @rdname as.matrix-methods
#' @export
setMethod(f="as.matrix", signature=signature(x = "SparseNeuroVec"), def=function(x,...) {
			  as(x, "matrix")
		  })

#' as.list
#'
#' convert SparseNeuroVec to list of \code{\linkS4class{DenseNeuroVol}}
#'
#' @rdname as.list-methods
#' @export
setMethod(f="as.list", signature=signature(x = "SparseNeuroVec"), def=function(x) {
			D4 <- dim(x)[4]
			lapply(1:D4, function(i) x[[i]])

})

#' show a \code{SparseNeuroVec}
#' @param object the object
#' @importFrom crayon bold blue green red yellow silver
#' @export
setMethod("show",
          signature=signature(object="SparseNeuroVec"),
          def=function(object) {
            # Get class name without package prefix
            class_name <- sub(".*:", "", class(object)[1])

            # Header
            cat("\n", crayon::bold(crayon::blue(class_name)), " ", sep = "")
            if (nchar(object@label) > 0) {
              cat(crayon::silver(paste0("'", object@label, "'")), "\n")
            } else {
              cat("\n")
            }

            # Spatial Info Block
            cat(crayon::bold("\n+= Spatial Info "), crayon::silver("---------------------------"), "\n", sep = "")
            dims <- dim(object)
            cat("| ", crayon::yellow("Dimensions"), "    : ", paste(dims[1:3], collapse = " x "), "\n", sep = "")
            cat("| ", crayon::yellow("Time Points"), "   : ", dims[4], "\n", sep = "")
            cat("| ", crayon::yellow("Spacing"), "       : ", paste(spacing(object)[1:3], collapse = " x "), "\n", sep = "")
            cat("| ", crayon::yellow("Origin"), "        : ", paste(round(origin(object)[1:3], 2), collapse = " x "), "\n", sep = "")

            # Sparse-specific Info Block
            card <- length(object@map@indices)
            cat(crayon::bold("\n+- Sparse Info  "), crayon::silver("----------------------------"), "\n", sep = "")
            cat("| ", crayon::yellow("Cardinality"), "   : ", card, "\n", sep = "")

            # Memory Usage Block
            mem_size <- object.size(object)
            size_str <- if (mem_size < 1024) {
              paste0(round(mem_size, 2), " B")
            } else if (mem_size < 1024^2) {
              paste0(round(mem_size/1024, 2), " KB")
            } else if (mem_size < 1024^3) {
              paste0(round(mem_size/1024^2, 2), " MB")
            } else {
              paste0(round(mem_size/1024^3, 2), " GB")
            }
            cat(crayon::bold("\n+= Memory Usage "), crayon::silver("--------------------------"), "\n", sep = "")
            cat("  ", crayon::yellow("Size"), "          : ", size_str, "\n", sep = "")
            cat("\n")
          })
