#' @include all_class.R
NULL
#' @include all_generic.R
NULL



#' NeuroVec: Neuroimaging Data Vector Class
#'
#' @description
#' The \code{NeuroVec} class represents a vectorized form of neuroimaging data, supporting both in-memory
#' and file-backed data modes. It provides efficient data storage and access methods and integrates with
#' the spatial reference system provided by \code{\linkS4class{NeuroSpace}}.
#'
#' @param data The image data. This can be:
#'   \itemize{
#'     \item A matrix (voxels x time points)
#'     \item A 4D array
#'     \item A list of \code{\linkS4class{NeuroVol}} objects
#'   }
#'   If a list of NeuroVol objects is provided, the geometric space (\code{NeuroSpace})
#'   will be inferred from the constituent volumes, which must all be identical.
#' @param space An optional \code{\linkS4class{NeuroSpace}} object defining the spatial
#'   properties of the image. Not required if \code{data} is a list of NeuroVol objects.
#' @param mask An optional logical array specifying which voxels to include. If provided,
#'   a SparseNeuroVec object will be created.
#' @param label A character string providing a label for the NeuroVec object. Default is an empty string.
#'
#' @return A concrete instance of the \code{\linkS4class{NeuroVec}} class:
#'   \itemize{
#'     \item If \code{mask} is provided: a \code{\linkS4class{SparseNeuroVec}} object
#'     \item Otherwise: a \code{\linkS4class{DenseNeuroVec}} object
#'   }
#'
#' @details
#' The function performs several operations:
#' \itemize{
#'   \item If \code{data} is a list of NeuroVol objects, it combines them into a single 4D array.
#'   \item It checks that the dimensions of \code{data} match the provided \code{space}.
#'   \item Depending on whether a \code{mask} is provided, it creates either a DenseNeuroVec
#'     or a SparseNeuroVec object.
#' }
#'
#' @examples
#' # Load an example 4D brain image
#' example_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
#' example_4d_image <- read_vec(example_file)
#'
#' # Create a DenseNeuroVec object
#' dense_vec <- NeuroVec(data = example_4d_image@.Data,
#'                       space = space(example_4d_image))
#' print(dense_vec)
#'
#' # Create a SparseNeuroVec object with a mask
#' mask <- array(runif(prod(dim(example_4d_image)[1:3])) > 0.5,
#'               dim = dim(example_4d_image)[1:3])
#' sparse_vec <- NeuroVec(data = example_4d_image@.Data,
#'                        space = space(example_4d_image),
#'                        mask = mask)
#' print(sparse_vec)
#'
#' @seealso \code{\linkS4class{NeuroSpace}} for spatial information,
#' \code{\link{sub_vector}} for subsetting routines, and
#' \code{\link{index_to_coord}} for coordinate conversion.
#' \code{\link{DenseNeuroVec-class}}, \code{\link{SparseNeuroVec-class}} for the
#' specific NeuroVec types.
#' \code{\link{NeuroVol-class}} for 3D volumetric data.
#'
#' @export
#' @importFrom methods new
#' @rdname NeuroVec-class
NeuroVec <- function(data, space = NULL, mask = NULL, label = "") {
  if (is.list(data)) {
    space <- space(data[[1]])
    space <- add_dim(space, length(data))
    data <- do.call(cbind, lapply(data, function(x) as.vector(x)))
  }


	if (is.null(mask)) {
	  if (prod(dim(space)) != length(data)) {
	    stop("dimensions of data argument do not match dimensions of space argument")
	  }
		DenseNeuroVec(data,space, label)
	} else {
		#SparseNeuroVec(data,space,mask,label)
	  SparseNeuroVec(data,space,mask)
	}

}

#' Create NeuroVec from list of NeuroVol objects
#'
#' @description
#' Factory function to create a NeuroVec object from a list of NeuroVol objects.
#' This is a convenience wrapper around the NeuroVec constructor that combines
#' multiple 3D volumes into a single 4D NeuroVec.
#'
#' @param vols A list of \code{\linkS4class{NeuroVol}} objects. All volumes must have
#'   identical spatial dimensions.
#' @param mask An optional logical array or \code{\linkS4class{LogicalNeuroVol}} object
#'   defining the subset of voxels to include. If provided, a SparseNeuroVec will be created.
#'
#' @return A \code{\linkS4class{NeuroVec}} object (either DenseNeuroVec or SparseNeuroVec
#'   depending on whether a mask is provided).
#'
#' @examples
#' # Create a simple NeuroVec from list of volumes
#' spc <- NeuroSpace(c(10, 10, 10))
#' vol1 <- NeuroVol(rnorm(10*10*10), spc)
#' vol2 <- NeuroVol(rnorm(10*10*10), spc)
#' vec <- vec_from_vols(list(vol1, vol2))
#' print(dim(vec))  # Should be c(10, 10, 10, 2)
#'
#' @export
#' @seealso \code{\link{NeuroVec}}, \code{\link{NeuroVol}}
vec_from_vols <- function(vols, mask = NULL) {
  if (!is.list(vols)) {
    stop("Input must be a list of NeuroVol objects")
  }
  
  if (length(vols) == 0) {
    stop("List of volumes cannot be empty")
  }
  
  # Check that all elements are NeuroVol objects
  if (!all(sapply(vols, function(x) inherits(x, "NeuroVol")))) {
    stop("All elements in the list must be NeuroVol objects")
  }
  
  # Check that all volumes have the same dimensions
  dims <- lapply(vols, dim)
  if (!all(sapply(dims, function(d) identical(d, dims[[1]])))) {
    stop("All volumes must have the same dimensions")
  }
  
  # Call the NeuroVec constructor with the list
  NeuroVec(vols, mask = mask)
}



#' DenseNeuroVec
#'
#' @description
#' This function constructs a DenseNeuroVec object, which represents a dense
#' four-dimensional brain image. It handles various input data formats and
#' ensures proper dimensionality.
#'
#' @param data The image data. This can be:
#'   \itemize{
#'     \item A 4-dimensional array
#'     \item A 2-dimensional matrix (either nvoxels x ntime-points or
#'           ntime-points x nvoxels)
#'     \item A vector (which will be reshaped to match the space dimensions)
#'   }
#' @param space A \code{\linkS4class{NeuroSpace}} object defining the spatial properties of the image.
#' @param label A character string providing a label for the DenseNeuroVec object. Default is an empty string.
#'
#' @return A concrete instance of the \code{\linkS4class{DenseNeuroVec}} class.
#'
#' @details
#' The function performs several operations based on the input data type:
#' \itemize{
#'   \item For matrix input: It determines the correct orientation (voxels x time or time x voxels)
#'     and reshapes accordingly. If necessary, it adds a 4th dimension to the space object.
#'   \item For vector input: It reshapes the data to match the dimensions specified in the space object.
#'   \item For array input: It ensures the dimensions match those specified in the space object.
#' }
#'
#' Note that the label parameter is currently not used in the object creation,
#' but is included for potential future use or consistency with other constructors.
#'
#' @examples
#' # Create a simple 4D brain image
#' dim <- c(64, 64, 32, 10)  # 64x64x32 volume with 10 time points
#' data <- array(rnorm(prod(dim)), dim)
#' space <- NeuroSpace(dim, spacing = c(3, 3, 4))
#'
#' # Create a DenseNeuroVec object
#' dense_vec <- DenseNeuroVec(data = data, space = space, label = "Example")
#' print(dense_vec)
#'
#' # Create from a matrix (voxels x time)
#' mat_data <- matrix(rnorm(prod(dim)), nrow = prod(dim[1:3]), ncol = dim[4])
#' dense_vec_mat <- DenseNeuroVec(data = mat_data, space = space)
#' print(dense_vec_mat)
#'
#' @seealso
#' \code{\link{NeuroVec-class}} for the parent class.
#' \code{\link{SparseNeuroVec-class}} for the sparse version of 4D brain images.
#' \code{\link{NeuroSpace-class}} for details on spatial properties.
#'
#' @export
#' @rdname DenseNeuroVec-class
DenseNeuroVec <- function(data, space, label="none") {
	if (is.matrix(data)) {
		splen <- prod(dim(space)[1:3])
		data <- if (nrow(data) == splen) {
			data
		} else if (ncol(data) == splen) {
			t(data)
		} else {
			stop("matrix dimensions do not match space dimensions")
		}

    if (length(dim(space)) == 3) {
      ## add 4th dim to space arg
      space <- add_dim(space, ncol(data))
    }

		dim(data) <- dim(space)
	} else if (is.vector(data) || length(dim(data)) == 1) {
	  data <- array(data, dim(space))
	}

	new("DenseNeuroVec", .Data=data, space=space, label=label)

}




#' Load image data from a NeuroVecSource object
#'
#' @description
#' This function loads the image data from a NeuroVecSource object, handling various
#' dimensionalities and applying any necessary transformations.
#'
#' @param x The NeuroVecSource object containing the image metadata and file information.
#'
#'
#' @details
#' This method performs the following steps:
#' 1. Validates the dimensionality of the metadata.
#' 2. Reads the image data using RNifti.
#' 3. Handles 5D arrays by dropping the 4th dimension if it has length 1.
#' 4. Applies slope scaling if present in the metadata.
#' 5. Constructs a NeuroSpace object with appropriate dimensions and spatial information.
#' 6. Creates and returns a DenseNeuroVec object, handling both 3D and 4D input arrays.
#'
#' @note This method currently only supports NIfTI file format through RNifti.
#'
#' @seealso \code{\link{NeuroVecSource}}, \code{\link{DenseNeuroVec}}, \code{\link{NeuroSpace}}
#'
#' @rdname load_data-methods
#' @return a DenseNeuroVec object
#' @export
setMethod(f="load_data", signature=c("NeuroVecSource"),
          def=function(x) {
            meta <- x@meta_info
            stopifnot(length(meta@dims) == 4)

            ind <- as.integer(x@indices)
            nels <- prod(meta@dims[1:3])

            is_gzip <- identical(meta@descriptor@data_encoding, "gzip") || endsWith(meta@data_file, ".gz")
            mmap_ok <- !is_gzip && identical(.Platform$endian, meta@endian)

            dat <- if (mmap_ok) {
              # Returns [time x voxels] for requested volumes
              read_mapped_vols(meta, ind)
            } else {
              # Stream volumes sequentially (works for gzipped inputs too).
              reader <- data_reader(meta, offset = 0)
              on.exit(close(reader), add = TRUE)

              pos <- split(seq_along(ind), ind)
              out <- matrix(0, nrow = length(ind), ncol = nels)
              max_vol <- max(ind)

              for (t in seq_len(max_vol)) {
                vol_dat <- read_elements(reader, nels)
                rows <- pos[[as.character(t)]]
                if (!is.null(rows)) {
                  out[rows, ] <- .apply_data_scaling(vol_dat, meta, index = t)
                }
              }
              out
            }

            # Apply scaling for mmap path (streaming path already scaled)
            if (mmap_ok) {
              dat <- .apply_data_scaling_matrix(dat, meta, indices = ind)
            }

            bspace <- NeuroSpace(c(meta@dims[1:3], length(ind)),
                                 meta@spacing,
                                 meta@origin,
                                 meta@spatial_axes,
                                 trans(meta))

            DenseNeuroVec(dat, bspace, label = meta@data_file)
          })






#' NeuroVecSource
#'
#' @description
#' This function constructs a NeuroVecSource object, which represents the source
#' of a four-dimensional brain image.
#'
#' @param file_name The name of the 4-dimensional image file.
#' @param indices An optional integer vector specifying the subset of volume indices to load.
#'   If not provided, all volumes will be loaded.
#' @param mask An optional logical array or \code{\linkS4class{NeuroVol}} object defining
#'   the subset of voxels to load. If provided, a SparseNeuroVecSource object will be created.
#'
#' @return An instance of the \code{\linkS4class{NeuroVecSource}} class.
#'
#' @details
#' If a \code{mask} is supplied, it should be a \code{\linkS4class{LogicalNeuroVol}} or
#' \code{\linkS4class{NeuroVol}} instance. If the latter, then the mask will be defined by
#' nonzero elements of the volume.
#'
#' @rdname NeuroVecSource
NeuroVecSource <- function(file_name, indices=NULL, mask=NULL) {
	if (!is.character(file_name)) {
	  cli::cli_abort("{.arg file_name} must be a character string.")
	}
	if (!file.exists(file_name)) {
	  cli::cli_abort("File {.path {file_name}} does not exist.")
	}

	meta_info <- read_header(file_name)

	if (!is.null(indices) && max(indices) > 1) {
	  if (length(dim(meta_info)) != 4) {
	    cli::cli_abort("Image must be 4-dimensional when {.arg indices} is specified, not {length(dim(meta_info))}D.")
	  }
	  if (max(indices) > dim(meta_info)[4]) {
	    cli::cli_abort("Max index {.val {max(indices)}} exceeds 4th dimension size {.val {dim(meta_info)[4]}}.")
	  }
	  if (min(indices) <= 0) {
	    cli::cli_abort("All indices must be > 0, but min is {.val {min(indices)}}.")
	  }
	}

  if (length(meta_info@dims) == 2) {
    stop(paste("cannot create NeuroVec with only two dimensions: ", paste(meta_info@dims, collapse=" ")))
  }

  if ( length(meta_info@dims) == 3) {
		indices <- 1
    meta_info@dims <- c(meta_info@dims,1)
	} else if (length(meta_info@dims) == 4 && is.null(indices)) {
		indices=seq(1, meta_info@dims[4])
	}

	if (is.null(mask)) {
		new("NeuroVecSource", meta_info=meta_info, indices=as.integer(indices))
	} else {
	  mask <- as.logical(mask)
	  if (any(is.na(mask))) {
	    cli::cli_abort("{.arg mask} must not contain NA values.")
	  }
		SparseNeuroVecSource(meta_info, as.integer(indices), mask)
	}

}




#' Get length of NeuroVec object
#'
#' @description
#' Returns the number of time points (4th dimension) in a NeuroVec object.
#' This represents the temporal dimension of the neuroimaging data.
#'
#' @param x A NeuroVec object
#'
#' @export
#' @rdname length-methods
setMethod("length", signature=c(x="NeuroVec"),
          def=function(x) {
            dim(x)[4]
          })



#' read_vol_list
#'
#' @description
#' This function loads a list of image volumes and returns a NeuroVec object.
#'
#' @param file_names A list of file names to load.
#' @param mask An optional mask defining the subset of voxels to load.
#'
#' @return An instance of the \code{\linkS4class{NeuroVec}} class.
#'
#' @export
#' @importFrom purrr map_lgl
read_vol_list <- function(file_names, mask=NULL) {
	stopifnot(all(map_lgl(file_names, file.exists)))
	meta_info <- lapply(file_names, read_header)

	dims <- do.call(rbind, lapply(meta_info, dim))
	if (!all(map_lgl(1:nrow(dims), function(i) all.equal(dims[1,], dims[i,])))) {
		stop("list of volumes must all have same dimensions")
	}

	if (!all(apply(dims, 1, length) == 3)) {
		stop("all volumes in list must have dim = 3")
	}

	nvols <- length(file_names)
	sourceList <- lapply(file_names, function(fname) {
		NeuroVolSource(fname, 1)
	})

	vols <- lapply(sourceList, load_data)
	if (is.null(mask)) {
		mat <- do.call(cbind, lapply(vols, function(v) as.vector(v@.Data)))
		dspace <- add_dim(space(vols[[1]]), length(vols))
		DenseNeuroVec(mat, dspace, label="")
	} else {
		mat <- do.call(cbind, lapply(vols, function(v) as.vector(v@.Data)))
		dspace <- add_dim(space(vols[[1]]), length(vols))
		if (is.vector(mask)) {
			## mask supplied as index vector, convert to logical
			M <- array(logical(prod(dim(dspace)[1:3])), dim(dspace)[1:3])
			M[mask] <- TRUE
			mask <- M
		} else {
			mask <- as.logical(mask)
		}


		#SparseNeuroVec(mat[mask,], dspace, mask=mask, label=map_chr(meta_info, function(m) m@label))
		SparseNeuroVec(mat[mask,], dspace, mask=mask)

	}
}

#' Drop a dimension
#' @export
#' @param x the object to drop a dimension from
#' @rdname drop-methods
#' @return An object of the same class as x with reduced dimensions or elements.
setMethod("drop", signature=(x="NeuroVec"),
          def=function(x) {
            if (dim(x)[4] == 1) {
              x[[1]]
            } else {
              x
            }
          })


# Coerce DenseNeuroVec to array.
setAs("DenseNeuroVec", "array", function(from) from@.Data)

# Coerce NeuroVec to array.
setAs("NeuroVec", "array", function(from) {
  vals <- from[]
  dim(vals) <- dim(from)
  vals
})

#' @rdname show-methods
#' @export
setMethod("show", "NeuroVecSource", function(object) {
  show_header("NeuroVecSource")
  show_field("Indices", paste(head(object@indices, 6), collapse = ", "),
             if (length(object@indices) > 6) " ..." else "")
  show(object@meta_info)
})




#' @rdname show-methods
#' @export
setMethod("show", "NeuroVec", function(object) {
  d <- dim(object)
  class_name <- sub(".*:", "", class(object)[1])
  show_header(class_name, format_mem(object))
  show_rule("Spatial")
  show_field("Dimensions", paste(d[1:3], collapse = " x "))
  show_field("Time Points", d[4])
  show_field("Spacing", paste(spacing(object)[1:3], collapse = " x "))
  show_field("Origin", paste(round(origin(object)[1:3], 2), collapse = ", "))
  show_field("Orientation", safe_axcodes(space(object)))
})



#' @export
#' @rdname concat-methods
setMethod(f="concat", signature=signature(x="NeuroVec", y="NeuroVol"),
		def=function(x,y,...) {
			.concat4D(x,y,...)
		})


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="NeuroVol", y="NeuroVec"),
  def=function(x,y,...) {
    .concat4D(x,y,...)
  })



#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="NeuroVec", y="NeuroVec"),
		def=function(x,y,...) {
			.concat4D(x,y,...)
		})



#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="ROIVec", y="ROIVec"),
          def=function(x,y,...) {
            ll <- list(x,y,...)

            cds <- map(ll, ~ coords(.))
            ident <- map_lgl(cds, ~ all(cds[[1]] == .))
            if (!all(ident)) {
              cli::cli_abort("concat.ROIVec: all {.cls ROIVec} arguments must have the same set of coordinates.")
            }

            dat <- do.call(rbind, map(ll, ~ .@.Data))
            vspace <- space(x)
            nd <- c(dim(vspace)[1:3], nrow(dat))


            nspace <-
              NeuroSpace(
                nd,
                origin = origin(x@space),
                spacing = spacing(x@space),
                axes = axes(x@space),
                trans = trans(x@space)
              )


            ROIVec(nspace, coords=coords(x), data=dat)
          })



#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="matrix"),
		def=function(x,i) {
			if (ncol(i) != 3) {
			  cli::cli_abort("Coordinate matrix {.arg i} must have 3 columns, not {ncol(i)}.")
			}

		  d4 <- dim(x)[4]
		  expanded <- i[rep(1:nrow(i), each=d4),]
		  expanded <- cbind(expanded, 1:d4)
	    vec <- x[expanded]
	    matrix(vec, d4, nrow(i))
    })


#' @rdname series-methods
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="matrix"),
          def=function(x,i) {
            mat <- series(x, i)
            ROIVec(space(x), coords=i, data=mat)

          })


#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="ROICoords"),
          def=function(x,i) {
            callGeneric(x, coords(i))
          })


#' @rdname series-methods
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="ROICoords"),
          def=function(x,i) {
            rvol <- series(x, coords(i))
            ROIVec(space(x), coords=coords(i), data=rvol)
          })


#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="LogicalNeuroVol"),
          def=function(x,i) {
            if (!isTRUE(all.equal(dim(x)[1:3], dim(i)[1:3]))) {
              cli::cli_abort("Spatial dimensions of {.arg x} ({.val {dim(x)[1:3]}}) and {.arg i} ({.val {dim(i)[1:3]}}) must match.")
            }
            idx <- which(i == TRUE)
            if (length(idx) == 0) {
              cli::cli_abort("{.arg i} mask contains no TRUE voxels.")
            }

            grid <- index_to_grid(i, idx)
            callGeneric(x, grid)

          })

#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="NeuroVol"),
          def=function(x,i) {
            callGeneric(x, as.logical(i))
          })

#' @rdname series-methods
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="LogicalNeuroVol"),
          def=function(x,i) {
            mat <- as.matrix(series(x, i))
            ROIVec(space(x), coords=index_to_grid(i, which(i == TRUE)), data=as.matrix(mat))
          })



#' @export
#' @param j second dimension index
#' @param k third dimension index
#' @param drop whether to drop dimension of length 1
#' @rdname series-methods
setMethod(f="series", signature(x="NeuroVec", i="integer"),
          def=function(x,i,j,k, drop=TRUE) {
            if (missing(j) && missing(k)) {
              nels <- prod(dim(x)[1:3])
              offsets <- seq(0, dim(x)[4]-1) * nels
              idx <- as.vector(outer(offsets, i, "+"))
              vals <- x[idx]
              ret <- matrix(vals, dim(x)[4], length(i))
              if (isTRUE(drop)) base::drop(ret) else ret
            } else {
              ## could be solved with expand.grid, no?
              if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
                cli::cli_abort("When providing {.arg j} and {.arg k}, {.arg i}, {.arg j}, and {.arg k} must each be length 1.")
              }
              ret <- x[i,j,k,]
              if (isTRUE(drop)) base::drop(ret) else ret
            }
          })


#' @rdname series-methods
#' @export
setMethod("series", signature(x="DenseNeuroVec", i="matrix"),
          def=function(x, i) {
            if (ncol(i) != 3) {
              cli::cli_abort("Coordinate matrix {.arg i} must have 3 columns, not {ncol(i)}.")
            }
            if (!is.numeric(i) || any(!is.finite(i))) {
              cli::cli_abort("Coordinate matrix {.arg i} must contain only finite numeric values.")
            }
            if (any(i != as.integer(i))) {
              cli::cli_abort("Coordinate matrix {.arg i} must contain whole-number voxel coordinates.")
            }
            d <- dim(x)
            validate_indices(d[1:3], list(i[,1], i[,2], i[,3]), c("i", "j", "k"))
            i <- matrix(as.integer(i), ncol = 3)
            # Direct linear indexing into .Data — avoids S4 dispatch overhead
            lin <- (i[,3] - 1L) * d[1] * d[2] + (i[,2] - 1L) * d[1] + i[,1]
            nt <- d[4]
            nels <- prod(d[1:3])
            out <- matrix(0, nt, nrow(i))
            for (t in seq_len(nt)) {
              out[t, ] <- x@.Data[lin + (t - 1L) * nels]
            }
            out
          })

#' @export
#' @param j second dimension index
#' @param k third dimension index
#' @param drop whether to drop dimension of length 1
#' @rdname series-methods
setMethod(f="series", signature=signature(x="DenseNeuroVec", i="integer"),
          def=function(x,i,j,k,drop=TRUE) {
            if (missing(j) && missing(k)) {
              g <- indexToGridCpp(i, dim(x)[1:3])
              ret <- callGeneric(x,g)
              if (isTRUE(drop)) base::drop(ret) else ret
            } else {
              if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
                cli::cli_abort("When providing {.arg j} and {.arg k}, {.arg i}, {.arg j}, and {.arg k} must each be length 1.")
              }
              ret <- x[i,j,k,]
              if (isTRUE(drop)) base::drop(ret) else ret
            }
          })



#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="numeric"),
		def=function(x, i, j, k, drop=TRUE) {
		  if (missing(j) && missing(k)) {
			  callGeneric(x,as.integer(i), drop=drop)
		  } else {
		    callGeneric(x,as.integer(i), as.integer(j), as.integer(k))
		  }
		})


#' @rdname series-methods
#' @param i first dimension index
#' @param j second dimension index
#' @param k third dimension index
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="numeric"),
          def=function(x, i, j, k, drop=TRUE) {
            mat <- if (missing(j) && missing(k)) {
              vdim <- dim(x)[1:3]
              vox <- arrayInd(i, vdim)
              callGeneric(x, vox)
            } else if (missing(i) || missing(j) || missing(k)) {
              stop("series_roi: must provide either 1D 'i' or 3D ('i', 'j', 'k') vector indices")
            } else {
              if (length(i) != 1 || length(j) != 1 || length(k) != 1) {
                cli::cli_abort("When providing {.arg j} and {.arg k}, {.arg i}, {.arg j}, and {.arg k} must each be length 1.")
              }
              ret <- x[i,j,k,]
              if (isTRUE(drop)) base::drop(ret) else ret

            }


          })







# Coerce NeuroVec to matrix.
setAs(from="NeuroVec", to="matrix",
      function(from) {
        dm <- dim(from)
        d123 <- prod(dm[1:3])
        d4 <- dm[4]
        vals <- from[]
        dim(vals) <- c(d123, d4)
        vals
      })


# Coerce DenseNeuroVec to matrix.
setAs(from="DenseNeuroVec", to="matrix",
		function(from) {
			data <- from@.Data
			dm <- dim(data)
			d123 <- prod(dm[1:3])
			d4 <- dm[4]

			dim(data) <- c(d123,d4)
			return(data)

		})


#' convert a \code{NeuroVec} to \code{list} of volumes.
#'
#' @rdname as.list-methods
#' @param x the object
#' @export
setMethod(f="as.list", signature=signature(x = "NeuroVec"), def=function(x) {
  map(1:length(x), ~ x[[.]])
})


#' convert a \code{NeuroVec} to a matrix
#'
#' @rdname as.matrix-methods
#' @param x the object
#' @param ... Additional arguments
#' @export
setMethod(f="as.matrix", signature=signature(x = "NeuroVec"), def=function(x) {
  dm <- dim(x)
  d123 <- prod(dm[1:3])
  d4 <- dm[4]
  vals <- as.vector(x@.Data)
  dim(vals) <- c(d123, d4)
  vals
})


# Coerce ROIVec to SparseNeuroVec.
setAs(from="ROIVec", to="SparseNeuroVec",
      function(from) {
        dat <- from@.Data
        mask <- array(0, dim(from@space)[1:3])
        mask[coords(from)] <- 1
        SparseNeuroVec(dat, from@space, mask=mask)
      })


#' @title Convert DenseNeuroVec to sparse representation using mask
#' @description This method converts a DenseNeuroVec object to a sparse representation using a given LogicalNeuroVol mask.
#' @param x A DenseNeuroVec object to convert to a sparse representation.
#' @param mask A LogicalNeuroVol object representing the mask to apply during conversion.
#' @return A SparseNeuroVec object resulting from the conversion.
#' @export
#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVec", mask="LogicalNeuroVol"),
          def=function(x, mask) {
            if (!all(dim(x)[1:3] == dim(mask))) {
              cli::cli_abort("Spatial dimensions of {.arg x} ({.val {dim(x)[1:3]}}) must match dimensions of {.arg mask} ({.val {dim(mask)}}).")
            }
            if (!all(spacing(x) == spacing(mask))) {
              cli::cli_abort("Spacing of {.arg x} and {.arg mask} must be identical.")
            }

            vdim <- dim(x)[1:3]
            dat <- as.matrix(x)[mask == TRUE,]
            bvec <- SparseNeuroVec(dat, space(x), mask)

          })


#' @title Convert DenseNeuroVec to sparse representation using a numeric mask
#' @description This method converts a DenseNeuroVec object to a sparse representation using a given numeric mask.
#' @param x A DenseNeuroVec object to convert to a sparse representation.
#' @param mask A numeric vector representing the mask to apply during conversion.
#' @return A SparseNeuroVec object resulting from the conversion.
#' @export
#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVec", mask="numeric"),
		def=function(x, mask) {
			vdim <- dim(x)[1:3]
			m <- array(0, vdim)
			m[mask] <- TRUE

			logivol <- LogicalNeuroVol(m, drop_dim(space(x)))

			dat <- as(x, "matrix")[mask,]

			bvec <- SparseNeuroVec(dat, space(x), logivol)

		})




#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec", signature=signature(x="NeuroHyperVec", file_name="character", format="missing", data_type="missing"),
          def=function(x, file_name) {
            write_nifti_hyper_vector(x, file_name)
          })

#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec", signature=signature(x="NeuroHyperVec", file_name="character", format="character", data_type="missing"),
          def=function(x, file_name, format, ...) {
            if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
              write_nifti_hyper_vector(x, file_name)
            } else {
              stop(paste("format ", format, "not supported for NeuroHyperVec."))
            }
          })

#' @export write_vec
#' @rdname write_vec-methods
#' @aliases write_vec,NeuroHyperVec,character,missing,character,ANY-method
setMethod(f="write_vec", signature=signature(x="NeuroHyperVec", file_name="character", format="missing", data_type="character"),
          def=function(x, file_name, data_type) {
            write_nifti_hyper_vector(x, file_name, data_type)
          })

#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="ROIVec", file_name="character", format="missing", data_type="missing"),
          def=function(x, file_name) {
            callGeneric(as(x, "SparseNeuroVec"), file_name)
          })


#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="missing"),
		def=function(x, file_name) {
			write_nifti_vector(x, file_name)
		})


#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="character", data_type="missing"),
		def=function(x, file_name, format, nbit=FALSE, compression=5, chunk_dim=c(10,10,10,dim(x)[4])) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, file_name)
			} else if (toupper(format) == "H5") {
			  if (!endsWith(file_name, ".h5")) {
			    file_name <- paste0(file_name, ".h5")
			  }

				h5obj <- to_nih5_vec(x, file_name, nbit=nbit, compression=compression, chunk_dim=chunk_dim)
				h5obj
			} else {
			  stop(paste("format ", format, "not supported."))
			}
		})


#' @export write_vec
#' @rdname write_vec-methods
#' @aliases write_vec,NeuroVec,character,missing,character,ANY-method
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="character"),
		def=function(x, file_name, data_type) {
			write_nifti_vector(x, file_name, data_type)

		})


#' @export
#' @rdname as_mmap
setMethod("as_mmap", signature(x = "NeuroVec"),
          function(x, file = NULL, data_type = "FLOAT", overwrite = FALSE, ...) {
            if (inherits(x, "MappedNeuroVec")) {
              return(x)
            }

            if (is.null(file)) {
              file <- tempfile(fileext = ".nii")
            } else if (file.exists(file) && !overwrite) {
              stop("as_mmap: file already exists; set overwrite = TRUE or choose a different path.")
            }

            if (endsWith(file, ".gz")) {
              stop("as_mmap: use an uncompressed .nii file for memory mapping.")
            }

            write_vec(x, file, data_type = data_type)
            MappedNeuroVec(file)
          })



#' @export
#' @rdname show-methods
setMethod("show", "DenseNeuroVec", function(object) {
  d <- dim(object)
  sp <- space(object)
  show_header("DenseNeuroVec", format_mem(object))
  show_rule("Spatial")
  show_field("Dimensions", paste(d[1:3], collapse = " x "),
             paste0(" (", d[4], " timepoints)"))
  show_field("Spacing", paste(spacing(sp)[1:3], collapse = " x "))
  show_field("Origin", paste(round(origin(sp)[1:3], 2), collapse = ", "))
  show_field("Orientation", safe_axcodes(sp))
  show_rule("Data")
  first_vol <- object[,,,1]
  show_field("Mean +/- SD", sprintf("%.3f +/- %.3f",
             mean(first_vol, na.rm = TRUE), sd(first_vol, na.rm = TRUE)),
             " (t=1)")
  if (!is.null(object@label) && nchar(object@label) > 0) {
    show_field("Label", object@label)
  }
})






#' @export
#' @rdname split_blocks-methods
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="integer"),
          def = function(x, indices,...) {
            if (length(indices) != dim(x)[4]) {
              cli::cli_abort("{.arg indices} length ({length(indices)}) must equal the 4th dimension of {.arg x} ({dim(x)[4]}).")
            }
            isplit <- split(1:length(indices), indices)

            f <- function(i) {
              sub_vector(x, isplit[[i]])
            }

            ret <- deflist::deflist(f, length(isplit), names=names(isplit))
            ret
          })

#' @export
#' @rdname split_blocks-methods
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="factor"),
          def = function(x, indices,...) {
            if (length(indices) != dim(x)[4]) {
              cli::cli_abort("{.arg indices} length ({length(indices)}) must equal the 4th dimension of {.arg x} ({dim(x)[4]}).")
            }
            ind <- as.integer(indices)
            ret <- callGeneric(x, ind)

            names(ret) <- levels(indices)
            ret
          })


#' read_vec
#'
#' @description
#' Loads a neuroimaging volume from one or more files, with support for various input formats
#' and memory management strategies.
#'
#' @param file_name The name(s) of the file(s) to load. If multiple files are specified,
#'   they are loaded and concatenated along the time dimension.
#' @param indices The indices of the sub-volumes to load (e.g. if the file is 4-dimensional).
#'   Only supported in "normal" mode.
#' @param mask A logical mask defining which spatial elements to load. Required for "bigvec" mode
#'   and optional for other modes.
#' @param mode The IO mode which is one of:
#'   * "normal": Standard in-memory loading
#'   * "mmap": Memory-mapped access (more memory efficient)
#'   * "bigvec": Optimized for large datasets with masking
#'   * "filebacked": File-backed storage with on-demand loading
#'
#' @details
#' This function supports multiple file formats:
#' * .nii: Standard NIfTI format
#' * .nii.gz: Compressed NIfTI (not supported in mmap mode)

#'
#' Memory management modes:
#' * "normal": Loads entire dataset into memory. Best for smaller datasets or when memory
#'   is not a constraint.
#' * "mmap": Memory-maps the file, providing efficient access for large files without
#'   loading entirely into memory. Not available for compressed files.
#' * "bigvec": Optimized for large datasets where only a subset of voxels are of interest.
#'   Requires a mask to specify which voxels to load.
#' * "filebacked": Similar to mmap but with more flexible caching strategies.
#'
#' @return An \code{\linkS4class{NeuroVec}} object representing the loaded volume(s).
#'
#' @examples
#'
#' # Load a single NIfTI file
#' img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#'
#' # Memory-mapped loading for large files
#' big_img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"), mode="mmap")
#'
#' # Load masked data for memory efficiency
#' mask <- as.logical(big_img[[1]])
#' masked_data <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"),
#'                mask=mask, mode="bigvec")
#'
#'
#' @export
#' @note
#' * Memory-mapping (.mmap mode) is not supported for gzipped files
#' * The bigvec mode requires a mask to be specified
#' * When loading multiple files, they must have compatible dimensions
read_vec  <- function(file_name, indices=NULL, mask=NULL, mode=c("normal", "mmap", "bigvec", "filebacked")) {
  mode <- match.arg(mode)

  get_source <- function(x, indices, mask) {
    if (endsWith(x, ".nii")) {
      NeuroVecSource(x, indices, mask)
    } else if (endsWith(x, ".nii.gz")) {
      NeuroVecSource(x, indices, mask)
    } #else if (endsWith(x, ".lv.h5")) {
      #if (!is.null(indices)) {
      #  warning("do not support indices argument for .lv.h5 files.")
      #}
      #if (!is.null(mask)) {
      #  warning("do not support mask argument for .lv.h5 files.")
      #}
      #
      #LatentNeuroVecSource(file_name)
    #} else if (endsWith(x, ".h5")) {
    #  if (!is.null(indices)) {
    #    warning("do not support indices argument for .h5 files.")
    #  }
    #  if (!is.null(mask)) {
    #    warning("do not support mask argument for .h5 files.")
    #  }
    #  H5NeuroVecSource(file_name=x)
    #}
  }

  vecs <- if (mode == "normal") {
	  lapply(file_name, function(fn) load_data(get_source(fn, indices, mask)))
	  #load_data(src)
  } else if (mode == "mmap") {
    if (!is.null(indices)) {
      stop("read_vec: memory mapped mode does not currently support volume 'indices'")
    }

    out <- vector(length(file_name), mode="list")

    for (i in seq_along(file_name)) {
      fname <- file_name[i]
      if (endsWith(fname, ".gz")) {
        if (!requireNamespace("R.utils")) {
          stop("must install 'R.utils' to enable temporary gzip decompression")
        }
        oname <- stringr::str_replace(fname, ".gz$", "")
        tdir <- tempdir()
        oname <- paste0(tdir, "/", basename(oname))
        R.utils::gunzip(fname, destname=oname, remove=FALSE)
        warning(paste("uncompressing gzipped file", fname, "to temporary file", oname, "for mem-mapping"))
        out[[i]] <- load_data(MappedNeuroVecSource(oname))
      } else {
        message("loading ", fname, " as mmaped file ")
        out[[i]] <- load_data(MappedNeuroVecSource(fname))
      }
    }
    out
  } else if (mode == "bigvec") {
    if (is.null(mask)) {
      stop("read_vec: 'bigvec' mode requires a mask")
    }

    out <- vector(length(file_name), mode="list")

    for (i  in seq_along(file_name)) {
      message("loading ", file_name[i], " as bigvec ")
      v <- load_data(NeuroVecSource(file_name[i], indices, mask))
      out[[i]] <- BigNeuroVec(v@data, space(v), mask)
    }

    out
  } else if (mode == "filebacked") {
    lapply(file_name, function(fn) FileBackedNeuroVec(fn))
  } else {
    stop()
  }

  if (length(vecs) == 1) {
    vecs[[1]]
  } else {
    do.call(NeuroVecSeq, vecs)
  }
}



 

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="integer"),
          def = function(x, clusters,...) {

            if (length(clusters) != prod(dim(x)[1:3])) {
              cli::cli_abort("{.arg clusters} length ({length(clusters)}) must equal the number of spatial voxels ({prod(dim(x)[1:3])}).")
            }
            keep <- which(clusters > 0 & !is.na(clusters))
            clusters <- clusters[keep]
            if (length(clusters) == 0) {
              cli::cli_abort("{.arg clusters} contains no positive, non-NA values.")
            }

            isplit <- split(keep, clusters)

            f <- function(i) {
              series_roi(x, isplit[[i]])
            }

            ret <- deflist::deflist(f, length(isplit), names=names(isplit))
            #ret <- deferred_list(out)
            #names(ret) <- names(isplit)
            ret
          })

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="numeric"),
          def = function(x, clusters,...) {
            callGeneric(x, as.integer(clusters), ...)
          })


#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="ClusteredNeuroVol"),
          def = function(x, clusters,...) {
            if (prod(dim(x)[1:3]) != length(clusters@mask)) {
              cli::cli_abort("Number of spatial voxels in {.arg x} ({prod(dim(x)[1:3])}) must equal length of {.arg clusters} mask ({length(clusters@mask)}).")
            }
            m <- which(clusters@mask > 0)
            clus <- rep(0, length(clusters@mask))
            clus[m] <- clusters@clusters
            split_clusters(x,clus)
          })

#' @rdname split_blocks-methods
#' @export
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="factor"),
          def = function(x, indices,...) {
            if (length(indices) != dim(x)[4]) {
              cli::cli_abort("{.arg indices} length ({length(indices)}) must equal the 4th dimension of {.arg x} ({dim(x)[4]}).")
            }
            ind <- as.integer(indices)
            ret <- callGeneric(x, ind)
            names(ret) <- levels(indices)
            ret
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="missing"),
          def = function(x) {
            ind <- 1:prod(dim(x)[1:3])
            time <- seq(1, dim(x)[4])
            lent = length(time)
            grid <- indexToGridCpp(ind, dim(x)[1:3])
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="DenseNeuroVec", subset="missing"),
          def = function(x) {
            ind <- 1:prod(dim(x)[1:3])
            time <- seq(1, dim(x)[4])
            lent <- length(time)
            grid <- indexToGridCpp(ind, dim(x)[1:3])
            f <- function(i) {
              imat <- cbind(do.call("rbind", rep(list(grid[i,]),lent)), time)
              x@.Data[imat]
            }
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="numeric"),
          def = function(x, subset) {
            ind <- subset
            if (max(ind) > prod(dim(x)[1:3])) {
              cli::cli_abort("Max index in {.arg subset} ({max(ind)}) exceeds number of spatial voxels ({prod(dim(x)[1:3])}).")
            }
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="logical"),
          def = function(x, subset) {
            if (length(subset) != prod(dim(x)[1:3])) {
              cli::cli_abort("{.arg subset} length ({length(subset)}) must equal the number of spatial voxels ({prod(dim(x)[1:3])}).")
            }
            ind <- which(subset)
            if (length(ind) == 0) {
              cli::cli_abort("{.arg subset} contains no TRUE values.")
            }
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vols-methods
setMethod(f="vols", signature=signature(x="NeuroVec", indices="numeric"),
          def = function(x, indices) {
            if (min(indices) <= 0 || max(indices) > dim(x)[4]) {
              cli::cli_abort("{.arg indices} must be in range [1, {dim(x)[4]}], got range [{min(indices)}, {max(indices)}].")
            }
            force(x)
            f <- function(i) x[[indices[i]]]
            #lis <- lapply(indices, function(i) function(i) x[[i]])
            #deferred_list(lis)
            deflist::deflist(f, length(indices))
          })

#' @export
#' @rdname vols-methods
setMethod(f="vols", signature=signature(x="NeuroVec", indices="missing"),
          def = function(x) {
            f <- function(i) x[[i]]
            #lis <- lapply(1:(dim(x)[4]), function(i) f)
            #deferred_list(lis)
            deflist::deflist(f, dim(x)[4])
          })

#' @rdname sub_vector-methods
#' @export
setMethod(f="sub_vector", signature=signature(x="NeuroVec", i="numeric"),
          def=function(x, i) {
            if (max(i) > dim(x)[4]) {
              cli::cli_abort("Max index {.val {max(i)}} exceeds 4th dimension size {.val {dim(x)[4]}}.")
            }
            xs <- space(x)
            dat <- x[,,,i]

            newdim <- c(dim(x)[1:3], length(i))
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVec(dat, bspace)
          })



#' @rdname sub_vector-methods
#' @export
setMethod(f="sub_vector", signature=signature(x="NeuroVecSeq", i="numeric"),
          def=function(x, i) {
            if (max(i) > dim(x)[4]) {
              cli::cli_abort("Max index {.val {max(i)}} exceeds 4th dimension size {.val {dim(x)[4]}}.")
            }
            lens <- sapply(x@vecs, function(v) dim(v)[4])
            offset <- c(0, cumsum(lens)) + 1

            vmap <- do.call(rbind, lapply(1:length(lens), function(j) {
              data.frame(i=j, offset=seq(offset[j], offset[j] + lens[j]-1), lind=1:lens[j])
            }))

            probe <- vmap[i,]
            smap <- split(probe$lind, probe$i)
            runs <- as.integer(names(smap))
            if (length(runs) == 0) {
              cli::cli_abort("No valid run segments found for the requested indices.")
            }

            svecs <- lapply(runs, function(rnum) {
              neuroim2::sub_vector(x@vecs[[rnum]], smap[[as.character(rnum)]])
            })

            out <- do.call(NeuroVecSeq, svecs)

            if (dim(out)[4] == 1) {
              ## simplify and return NeuroVec. Could further simplify by return NeuroVec if length(vecs) == 1
              out@vecs[[1]]
            } else {
              out
            }

          })



#' [[
#'
#' @description
#' This function extracts a single volume from a NeuroVec object.
#'
#' @param x The NeuroVec object.
#' @param i The volume index to extract.
#'
#' @return a DenseNeuroVol object
#' @rdname extract-methods
#' @export
setMethod(f="[[", signature=signature(x="NeuroVec", i="numeric"),
          def = function(x, i) {
            ## or ... drop(sub_vector(x,i))
            if (length(i) != 1) {
              cli::cli_abort("{.arg i} must be a single index, not length {length(i)}.")
            }
            xs <- space(x)
            # Use drop=FALSE to prevent dropping dimensions
            dat <- x[,,,i, drop=FALSE]
            newdim <- dim(x)[1:3]
            dat <- array(dat, dim=newdim)
            bspace <- NeuroSpace(newdim, spacing=spacing(xs),
                                 origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVol(dat, bspace)
          })


#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="ROIVec", file_name="character", format="missing", data_type="missing"),
          def=function(x, file_name) {
            callGeneric(as(x, "SparseNeuroVec"), file_name)
          })


#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="missing"),
		def=function(x, file_name) {
			write_nifti_vector(x, file_name)
		})


#' @export
#' @rdname write_vec-methods
#' @param nbit set nbit compression
#' @param compression compression level 1 to 9
#' @param chunk_dim the dimensions of each chunk
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="character", data_type="missing"),
		def=function(x, file_name, format, nbit=FALSE, compression=5, chunk_dim=c(10,10,10,dim(x)[4])) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, file_name)
			} else {
			  stop(paste("format ", format, "not supported."))
			}
		})


#' @export write_vec
#' @rdname write_vec-methods
#' @aliases write_vec,NeuroVec,character,missing,character,ANY-method
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="character"),
		def=function(x, file_name, data_type) {
			write_nifti_vector(x, file_name, data_type)

		})



#' @export
#' @rdname show-methods
setMethod("show", "NeuroVecSeq", function(object) {
  n <- length(object@vecs)
  show_header("NeuroVecSeq", paste(n, "vectors"))
  show_rule("Sequence")
  for (i in seq_len(min(n, 5))) {
    v <- object@vecs[[i]]
    d <- dim(v)
    show_field(paste0("  [", i, "]"), paste0(class(v)[1], " ",
               paste(d[1:3], collapse="x"), " x ", d[4], "t"))
  }
  if (n > 5) {
    cat("  ... and", n - 5, "more\n")
  }
})


#' @export
#' @rdname as.matrix-methods
setMethod("as.matrix", "DenseNeuroVec",
  function(x) {
    d <- dim(x)
    matrix(as.array(x@.Data), nrow = prod(d[1:3]), ncol = d[4])
  }
)

#' @rdname mask-methods
#' @export
setMethod("mask", "DenseNeuroVec",
          function(x) {
            # Extract 3D spatial dimensions for the mask
            spatial_dims <- dim(x)[1:3]
            LogicalNeuroVol(array(TRUE, spatial_dims), 
                           NeuroSpace(spatial_dims, 
                                     spacing(x)[1:3],
                                     origin(x)[1:3]))
          })
