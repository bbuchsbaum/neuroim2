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
#' @importFrom assertthat assert_that
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
	} else if (is.vector(data) || length(dim(data) == 1)) {
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
#' @return A DenseNeuroVec object containing the loaded image data and associated spatial information.
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
#' @examples
#' \dontrun{
#' # Assuming 'source' is a valid NeuroVecSource object
#' vec_data <- load_data(source)
#' }
#'
#' @rdname load_data-methods
#' @export
setMethod(f="load_data", signature=c("NeuroVecSource"),
          def=function(x) {
            #browser()
            meta <- x@meta_info

            stopifnot(length(meta@dims) == 4)

            nels <- prod(meta@dims[1:4])
            ind <- x@indices

            ## use RNifti, fails to work with other formats, though...
            arr <- RNifti::readNifti(meta@data_file)

            if (length(dim(arr)) == 5 && dim(arr)[4] == 1) {
              ## if 4th dimension is of length 1, drop it
              arr <- drop(arr)
            }


            ## bit of a hack to deal with scale factors
            if (.hasSlot(meta, "slope")) {

              if (meta@slope != 0) {
                arr <- arr * meta@slope
              }
            }

            bspace <- NeuroSpace(c(meta@dims[1:3], length(ind)),meta@spacing, meta@origin,
                                 meta@spatial_axes, trans(meta))

            if (length(dim(arr)) == 3) {
              dim(arr) <- c(dim(arr),1)
              DenseNeuroVec(unclass(arr), bspace, label=meta@data_file)
            } else if (length(dim(arr)) == 4) {
              DenseNeuroVec(arr[,,,ind,drop=FALSE], bspace, label=meta@data_file)
            } else {
              stop("NeuroVecSource::load_data: array dimension must be equal to 3 or 4.")
            }
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
#' @importFrom assertthat assert_that
NeuroVecSource <- function(file_name, indices=NULL, mask=NULL) {
	assert_that(is.character(file_name))
	assert_that(file.exists(file_name))


	meta_info <- read_header(file_name)

	if (!is.null(indices) && max(indices) > 1) {
	  assert_that(length(dim(meta_info)) == 4)
	  assert_that(max(indices) <= dim(meta_info)[4])
	  assert_that(min(indices) > 0)
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
	  assert_that(!any(is.na(mask)))
		SparseNeuroVecSource(meta_info, as.integer(indices), mask)
	}

}




#' Get length of NeuroVec object
#'
#' @description Returns the number of time points (4th dimension) in a NeuroVec object
#'
#' @param x A NeuroVec object
#' @return Integer length (number of time points)
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
		mat <- do.call(cbind, vols)
		dspace <- add_dim(space(vols[[1]]), length(vols))
		DenseNeuroVec(mat, dspace, label=map_chr(meta_info, function(m) m@label))
	} else {
		mat <- do.call(cbind, vols)
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

#' drop
#'
#' @description
#' This function drops the last dimension of a NeuroVec object if it is of length 1.
#'
#' @param x The NeuroVec object.
#'
#' @return If the last dimension of the NeuroVec object is of length 1, a DenseNeuroVol
#'   object is returned. Otherwise, the original NeuroVec object is returned.
#'
#' @export
setMethod("drop", signature=(x="NeuroVec"),
          def=function(x) {
            if (dim(x)[4] == 1) {
              idx <- seq(1, prod(dim(x)[1:3]))
              series(x, idx)
            } else {
              x
            }
          })


#' @export
setAs("DenseNeuroVec", "array", function(from) from@.Data)

#' @export
setAs("NeuroVec", "array", function(from) {
  vals <- from[]
  dim(vals) <- dim(from)
  vals
})

#' show a NeuroVecSource object
#'
#' @description
#' This function prints a summary of a NeuroVecSource object.
#'
#' @param object The NeuroVecSource object to display.
#'
#' @export
setMethod(f="show",
		signature=signature(object="NeuroVecSource"),
		def=function(object) {
			cat("an instance of class",  class(object), "\n\n")
			cat("   indices: ", object@indices, "\n\n")
			cat("   meta_info: \n")
			show(object@meta_info)
			cat("\n\n")

		})




#' show a NeuroVec object
#'
#' @description
#' This function prints a summary of a NeuroVec object.
#'
#' @param object The NeuroVec object to display.
#'
#' @export
setMethod(f="show", signature=signature("NeuroVec"),
          def=function(object) {
            # Get class name without package prefix
            class_name <- sub(".*:", "", class(object)[1])

            # Header
            cat("\n", crayon::bold(crayon::blue(class_name)), " ", sep="")
            if (nchar(object@label) > 0) {
              cat(crayon::silver(paste0("'", object@label, "'")), "\n", sep="")
            } else {
              cat("\n")
            }

            # Spatial Info
            cat(crayon::bold("\n- Spatial Info"), crayon::silver(" ---------------------------"), "\n", sep="")
            cat("| ", crayon::yellow("Dimensions"), "    : ", paste(dim(object)[1:3], collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Time Points"), "   : ", dim(object)[4], "\n", sep="")
            cat("| ", crayon::yellow("Spacing"), "       : ", paste(spacing(object)[1:3], collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Origin"), "        : ", paste(round(origin(object)[1:3], 2), collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Orientation"), "   : ", paste(object@space@axes@i@axis, object@space@axes@j@axis, object@space@axes@k@axis), "\n", sep="")

            # Memory Info
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

            cat(crayon::bold("\n- Memory Usage"), crayon::silver(" --------------------------"), "\n", sep="")
            cat("  ", crayon::yellow("Size"), "          : ", size_str, "\n", sep="")
            cat("\n")
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
            assert_that(all(ident), msg=paste("concat.ROIVec: ", "all 'ROIVec' arguments must have the same set of coordinates"))

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
			assertthat::assert_that(ncol(i) == 3)

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
            assertthat::assert_that(all.equal(dim(x)[1:3], dim(i)[1:3]))
            idx <- which(i == TRUE)
            assertthat::assert_that(length(idx) > 0)

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
              idx <- map(i, ~ . + offsets) %>% flatten_dbl()
              vals <- x[idx]
              ret <- matrix(vals, dim(x)[4], length(i))
              if (drop) drop(ret) else ret
            } else {
              ## could be solved with expand.grid, no?
              assert_that(length(i) == 1 && length(j) == 1 && length(k) ==1)
              ret <- x[i,j,k,]
              if (drop) drop(ret) else ret
            }
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
              if (drop) drop(ret) else ret
            } else {
              assert_that(length(i) == 1 && length(j) == 1 && length(k) ==1)
              ret <- x[i,j,k,]
              if (drop) drop(ret) else ret
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
          def=function(x, i, j, k) {
            mat <- if (missing(j) && missing(k)) {
              vdim <- dim(x)[1:3]
              vox <- arrayInd(i, vdim)
              callGeneric(x, vox)
            } else if (missing(i) || missing(j) || missing(k)) {
              stop("series_roi: must provide either 1D 'i' or 3D ('i', 'j', 'k') vector indices")
            } else {
              assert_that(length(i) == 1 && length(j) == 1 && length(k) ==1)
              ret <- x[i,j,k,]
              if (drop) drop(ret) else ret

            }


          })







#' @export
setAs(from="NeuroVec", to="matrix",
      function(from) {
        dm <- dim(from)
        d123 <- prod(dm[1:3])
        d4 <- dm[4]
        vals <- from[]
        dim(vals) <- c(d123, d4)
        vals
      })


#' @export
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


#' @export
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
            assert_that(all(dim(x)[1:3] == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))

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
#' @rdname show-methods
#' @aliases show,DenseNeuroVec-method
setMethod("show", "DenseNeuroVec",
          def=function(object) {
            # Get class name without package prefix
            class_name <- sub(".*:", "", class(object)[1])

            # Header with class name and memory info
            total_elements <- prod(dim(object))
            mem_size <- format(object.size(object) / 1024^2, digits=2)

            cat("\n", crayon::bold(crayon::blue("DenseNeuroVec")), " ",
                crayon::silver(paste0("(", mem_size, " MB)")), "\n", sep="")

            # Spatial information
            cat(crayon::bold("\n- Spatial Info"), crayon::silver(" ---------------------------"), "\n", sep="")
            cat("| ", crayon::yellow("Dimensions"), "    : ",
                paste(dim(object)[1:3], collapse=" x "),
                crayon::silver(" ("), dim(object)[4], " timepoints", crayon::silver(")"), "\n", sep="")
            cat("| ", crayon::yellow("Total Voxels"), "  : ",
                format(prod(dim(object)[1:3]), big.mark=","), "\n", sep="")
            cat("| ", crayon::yellow("Spacing"), "       : ",
                paste(object@space@spacing[1:3], collapse=" x "), "\n", sep="")

            # Data properties
            cat(crayon::bold("\n- Properties"), crayon::silver(" ---------------------------"), "\n", sep="")
            cat("| ", crayon::yellow("Origin"), "        : ",
                paste(round(object@space@origin[1:3], 2), collapse=" x "), "\n", sep="")
            cat("| ", crayon::yellow("Orientation"), "   : ",
                paste(object@space@axes@i@axis, object@space@axes@j@axis, object@space@axes@k@axis), "\n", sep="")

            # Summary statistics
            cat(crayon::bold("\n- Statistics"), crayon::silver(" ---------------------------"), "\n", sep="")

            # Calculate stats efficiently for the first timepoint
            first_vol <- object[,,,1]
            cat("    ", crayon::silver("Mean +/- SD"), "    : ",
                round(mean(first_vol, na.rm=TRUE), 3), " +/- ",
                round(sd(first_vol, na.rm=TRUE), 3), "\n", sep="")

            if(!is.null(object@label) && nchar(object@label) > 0) {
              cat("\n", crayon::italic(paste("Label:", object@label)), "\n", sep="")
            }

            cat("\n")
          })






#' @export
#' @rdname split_blocks-methods
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="integer"),
          def = function(x, indices,...) {
            assert_that(length(indices) == dim(x)[4])
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
            assert_that(length(indices) == dim(x)[4])
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
#' \dontrun{
#' # Load a single NIfTI file
#' img <- read_vec("subject01.nii")
#'
#' # Load multiple volumes and concatenate
#' imgs <- read_vec(c("run1.nii", "run2.nii"))
#'
#' # Memory-mapped loading for large files
#' big_img <- read_vec("large_volume.nii", mode="mmap")
#'
#' # Load masked data for memory efficiency
#' mask <- read_vol("brain_mask.nii")
#' masked_data <- read_vec("functional.nii", mask=mask, mode="bigvec")
#' }
#'
#' @export
#' @note
#' * Memory-mapping (.mmap mode) is not supported for gzipped files
#' * For .lv.h5 and .h5 files, the indices and mask parameters are ignored
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
#' @rdname as.matrix-methods
setMethod(f="as.matrix", signature=signature(x = "NeuroVec"), def=function(x) {
  dm <- dim(x)
  d123 <- prod(dm[1:3])
  d4 <- dm[4]
  vals <- as.vector(x@.Data)
  dim(vals) <- c(d123, d4)
  vals
})

#' @rdname series-methods
#' @param i first dimension index
#' @param j second dimension index
#' @param k third dimension index
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="numeric"),
          def=function(x, i, j, k) {
            mat <- if (missing(j) && missing(k)) {
              vdim <- dim(x)[1:3]
              vox <- arrayInd(i, vdim)
              callGeneric(x, vox)
            } else if (missing(i) || missing(j) || missing(k)) {
              stop("series_roi: must provide either 1D 'i' or 3D ('i', 'j', 'k') vector indices")
            } else {
              assert_that(length(i) == 1 && length(j) == 1 && length(k) ==1)
              ret <- x[i,j,k,]
              if (drop) drop(ret) else ret
            }
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
setMethod("series_roi", signature(x="NeuroVec", i="LogicalNeuroVol"),
          def=function(x,i) {
            mat <- as.matrix(series(x, i))
            ROIVec(space(x), coords=index_to_grid(i, which(i == TRUE)), data=as.matrix(mat))
          })

#' @export
#' @rdname split_clusters-methods
setMethod(f="split_clusters", signature=signature(x="NeuroVec", clusters="integer"),
          def = function(x, clusters,...) {

            assert_that(length(clusters) == prod(dim(x)[1:3]))
            keep <- which(clusters > 0 & !is.na(clusters))
            clusters <- clusters[keep]
            assert_that(length(clusters) > 0)

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
            assert_that(prod(dim(x)[1:3]) == length(clusters@mask))
            m <- which(clusters@mask > 0)
            clus <- rep(0, length(clusters@mask))
            clus[m] <- clusters@clusters
            split_clusters(x,clus)
          })

#' @rdname split_blocks-methods
#' @export
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="factor"),
          def = function(x, indices,...) {
            assert_that(length(indices) == dim(x)[4])
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
            assert_that(max(ind) <= prod(dim(x)[1:3]))
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="logical"),
          def = function(x, subset) {
            assert_that(length(subset) == prod(dim(x)[1:3]))
            ind <- which(subset)
            assert_that(length(ind) > 0)
            f <- function(i) series(x, ind[i])
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vols-methods
setMethod(f="vols", signature=signature(x="NeuroVec", indices="numeric"),
          def = function(x, indices) {
            assert_that(min(indices) > 0 && max(indices) <= dim(x)[4])
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
            assertthat::assert_that(max(i) <= dim(x)[4])
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
            assertthat::assert_that(max(i) <= dim(x)[4])
            lens <- sapply(x@vecs, function(v) dim(v)[4])
            offset <- c(0, cumsum(lens)) + 1

            vmap <- do.call(rbind, lapply(1:length(lens), function(j) {
              data.frame(i=j, offset=seq(offset[j], offset[j] + lens[j]-1), lind=1:lens[j])
            }))

            probe <- vmap[i,]
            smap <- split(probe$lind, probe$i)
            runs <- as.integer(names(smap))
            assertthat::assert_that(length(runs) > 0)

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
#' @return A DenseNeuroVol object representing the extracted volume.
#'
#' @export
setMethod(f="[[", signature=signature(x="NeuroVec", i="numeric"),
          def = function(x, i) {
            ## or ... drop(sub_vector(x,i))
            assert_that(length(i) == 1)
            xs <- space(x)
            dat <- x[,,,i]
            newdim <- dim(x)[1:3]
            bspace <- NeuroSpace(newdim, spacing=spacing(xs),
            origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVol(dat, bspace)
          })


#' drop
#'
#' @description
#' This function drops the last dimension of a NeuroVec object if it is of length 1.
#'
#' @param x The NeuroVec object.
#'
#' @return If the last dimension of the NeuroVec object is of length 1, a DenseNeuroVol
#'   object is returned. Otherwise, the original NeuroVec object is returned.
#'
#' @export
setMethod("drop", signature(x="NeuroVec"),
          def=function(x) {
            if (dim(x)[4] == 1) {
              idx <- seq(1, prod(dim(x)[1:3]))
              vals <- x[idx]
              sp <- drop_dim(space(x))
              DenseNeuroVol(array(vals, dim(sp)), sp)
            }
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
            assert_that(all(dim(x)[1:3] == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))

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
setMethod("show", "NeuroVecSeq",
          def=function(object) {
            cat("\n", crayon::bold(crayon::blue("NeuroVecSeq")), " ",
                crayon::silver(paste0("(", length(object@vecs), " vectors)")), "\n", sep="")

            cat(crayon::bold("\n- Sequence Info"), crayon::silver(" ---------------------------"), "\n", sep="")
            cat("  ", crayon::yellow("Length"), "        : ", length(object@vecs), "\n", sep="")
            cat("  ", crayon::yellow("Total Time"), "    : ", sum(object@lens), " points\n", sep="")

            sp <- space(object@vecs[[1]])
            cat(crayon::bold("\n- Spatial Info"), crayon::silver(" ---------------------------"), "\n", sep="")
            cat("  ", crayon::yellow("Dimensions"), "    : ", paste(dim(object@vecs[[1]])[1:3], collapse=" x "), "\n", sep="")
            cat("  ", crayon::yellow("Spacing"), "       : ", paste(sp@spacing[1:3], collapse=" x "), "\n", sep="")
            cat("  ", crayon::yellow("Origin"), "        : ", paste(round(sp@origin[1:3], 2), collapse=" x "), "\n", sep="")
            cat("  ", crayon::yellow("Orientation"), "   : ", paste(sp@axes@i@axis, sp@axes@j@axis, sp@axes@k@axis), "\n", sep="")

            cat(crayon::bold("\n- Vector Details"), crayon::silver(" --------------------------"), "\n", sep="")
            for (i in seq_along(object@vecs)) {
              v <- object@vecs[[i]]
              vclass <- sub(".*:", "", class(v)[1])
              cat("  ", crayon::green(paste0(i, ".")), " ",
                  crayon::cyan(vclass), " ",
                  crayon::silver(paste0("(", dim(v)[4], " timepoints)")),
                  "\n", sep="")
            }
            cat("\n")
          })


#' @export
#' @rdname as.matrix-methods
setMethod("as.matrix", "DenseNeuroVec",
  function(x) {
    d <- dim(x)
    matrix(as.array(x@.Data), nrow = prod(d[1:3]), ncol = d[4])
  }
)

