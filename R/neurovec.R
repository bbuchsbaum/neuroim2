#' @include all_class.R
NULL
#' @include all_generic.R
NULL



#' NeuroVec
#'
#' @description
#' This function constructs a NeuroVec object, which represents a four-dimensional 
#' brain image. It can create either a DenseNeuroVec or SparseNeuroVec object 
#' depending on the input parameters.
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
#' @seealso 
#' \code{\link{DenseNeuroVec-class}}, \code{\link{SparseNeuroVec-class}} for the 
#' specific NeuroVec types.
#' \code{\link{NeuroVol-class}} for 3D volumetric data.
#' \code{\link{NeuroSpace-class}} for details on spatial properties.
#'
#' @export
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
#' space <- NeuroSpace(dim, spacing = c(3, 3, 4, 1))
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
DenseNeuroVec <- function(data, space, label="") {

	if (is.matrix(data)) {
		splen <- prod(dim(space)[1:3])
		data <- if (ncol(data) == splen) {
			t(data)
		} else if (nrow(data) == splen) {
			data
		}

    if (length(dim(space)) == 3) {
      ## add 4th dim to space arg
      space <- add_dim(space, ncol(data))
    }

		dim(data) <- dim(space)
	} else if (is.vector(data) || length(dim(data) == 1)) {
	  data <- array(data, dim(space))
	}


	new("DenseNeuroVec", .Data=data,  space=space)

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
#' @noRd
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
              DenseNeuroVec(unclass(arr), bspace, x)
            } else if (length(dim(arr)) == 4) {
              DenseNeuroVec(arr[,,,ind,drop=FALSE], bspace, x)
            } else {
              stop("NeuroVecSource::load_data: array dimension must be equal to 3 or 4.")
            }
          })


#' Load data from an H5NeuroVecSource object
#'
#' @description
#' This function loads the image data from an H5NeuroVecSource object, which represents
#' data stored in an HDF5 file format.
#'
#' @param x The H5NeuroVecSource object containing the file name of the HDF5 data.
#'
#' @return An H5NeuroVec object containing the loaded image data.
#'
#' @details
#' This method simply calls the H5NeuroVec constructor with the file name stored in the source object.
#'
#' @seealso \code{\link{H5NeuroVecSource}}, \code{\link{H5NeuroVec}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'source' is a valid H5NeuroVecSource object
#' h5_data <- load_data(source)
#' }
#'
#' @noRd
setMethod(f="load_data", signature=c("H5NeuroVecSource"),
          def=function(x) {
            H5NeuroVec(x@file_name)
          })

#' Load data from a LatentNeuroVecSource object
#'
#' @description
#' This function loads and constructs a LatentNeuroVec object from a LatentNeuroVecSource,
#' which represents latent space neuroimaging data stored in an HDF5 file.
#'
#' @param x The LatentNeuroVecSource object containing the file name of the HDF5 data.
#'
#' @return A LatentNeuroVec object containing the loaded latent space representation of the neuroimaging data.
#'
#' @details
#' This method performs the following steps:
#' 1. Opens the HDF5 file specified in the source object.
#' 2. Reads the basis, loadings, offset, and indices data from the file.
#' 3. Constructs a NeuroSpace object from the spatial information in the file.
#' 4. Creates a mask from the indices data.
#' 5. Closes the HDF5 file.
#' 6. Constructs and returns a LatentNeuroVec object with the loaded data.
#'
#' @note This method requires the hdf5r package for reading HDF5 files.
#'
#' @seealso \code{\link{LatentNeuroVecSource}}, \code{\link{LatentNeuroVec}}, \code{\link{NeuroSpace}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'source' is a valid LatentNeuroVecSource object
#' latent_data <- load_data(source)
#' }
#'
#' @importFrom hdf5r H5File
#' @noRd
setMethod(f="load_data", signature=c("LatentNeuroVecSource"),
          def=function(x) {
            h5obj <- hdf5r::H5File$new(x@file_name)
            basis <- h5obj[["data/basis"]][,]
            loadings <- h5obj[["data/loadings"]][,]
            offset <- h5obj[["data/offset"]][]
            indices <- h5obj[["data/indices"]][]
            sp <- NeuroSpace(dim=h5obj[["space/dim"]][],
                             spacing=h5obj[["space/spacing"]][],
                             origin=h5obj[["space/origin"]][],
                             trans=h5obj[["space/trans"]][,])

            mask <- NeuroVol(array(0, dim(sp)[1:3]), drop_dim(sp))
            mask[indices] <- 1
            mask <- as.logical(mask)
            h5obj$close_all()
            LatentNeuroVec(basis, loadings, space=sp, mask=mask, offset=offset)
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




#' @keywords internal
#' @noRd
H5NeuroVecSource <- function(file_name) {
  new("H5NeuroVecSource", file_name=file_name)
}

#' H5NeuroVec
#'
#' @description
#' This function constructs an H5NeuroVec object, which represents a four-dimensional 
#' brain image stored in the neuroim2 HDF5 format.
#'
#' @param file_name The name of the 4-dimensional image file in neuroim2 HDF5 format.
#'
#' @return An instance of the \code{\linkS4class{H5NeuroVec}} class.
#'
#' @rdname NeuroVec
#' @importFrom assertthat assert_that
#' @export
H5NeuroVec <- function(file_name) {
  assert_that(is.character(file_name))
  assert_that(file.exists(file_name))


  h5obj <- hdf5r::H5File$new(file_name)

  rtype <- try(hdf5r::h5attr(h5obj, which="rtype"))
  if (! (rtype == "DenseNeuroVec")) {
    stop("invalid h5 file: ", file_name)
  }

  if (length(h5obj[["space/dim"]][]) != 4) {
    stop(paste("cannot H5NeuroVec: must have 4 dimensions: ", paste(h5obj[["space/dim"]][], collapse=" ")))
  }


  sp <- NeuroSpace(dim=h5obj[["space/dim"]][], origin=h5obj[["space/origin"]][],
                   trans=h5obj[["space/trans"]][,])

  new("H5NeuroVec", space=sp, obj=h5obj)

}


#' Get the length of a NeuroVec object
#'
#' @description
#' This function returns the number of volumes in a NeuroVec object.
#'
#' @param x The NeuroVec object.
#'
#' @return The number of volumes in the NeuroVec object.
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
setMethod("drop", signature(x="NeuroVec"),
          def=function(x) {
            if (dim(x)[4] == 1) {
              idx <- seq(1, prod(dim(x)[1:3]))
              vals <- x[idx]
              sp <- drop_dim(space(x))
              DenseNeuroVol(array(vals, dim(sp)), sp)
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
            sp <- space(object)
            cat(class(object), "\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(sp@spacing[1:(length(sp@spacing)-1)], " X ", collapse=" "),
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(sp@origin[1:(length(sp@origin)-1)], " X ", collapse=" "),
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")
            cat("  Coordinate Transform :", zapsmall(sp@trans), "\n")

          }
)






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

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="missing"),
          def = function(x) {
            ind <- 1:prod(dim(x)[1:3])
            time <- seq(1, dim(x)[4])
            lent = length(time)
            grid <- indexToGridCpp(ind, dim(x)[1:3])
            ##vox <- index_to_grid(x, ind)
            ##f <- function(i) series(x, vox[i,1], vox[i,2], vox[i,3])
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

            #f <- function(i) series(x, ind[i])
            f <- function(i) {
              imat <- cbind(do.call("rbind", rep(list(grid[i,]),lent)), time)
              x@.Data[imat]
            }
            #lis <- map(ind, function(i) f)
            #deferred_list(lis)
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="numeric"),
          def = function(x, subset) {
            ind <- subset
            assert_that(max(ind) <= prod(dim(x)[1:3]))

            ## 03-04-2023 index_to_grid is not necessary here
            ## vox <- index_to_grid(x, ind)
            ###############################################

            #f <- function(i) series(x, vox[i,1], vox[i,2], vox[i,3])
            f <- function(i) series(x, ind[i])

            #lis <- lapply(seq_along(ind), function(i) f)
            #deferred_list(lis)
            deflist::deflist(f, length(ind))
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="logical"),
          def = function(x, subset) {
            assert_that(length(subset) == prod(dim(x)[1:3]))
            ind <- which(subset)
            assert_that(length(ind) > 0)
            #vox <- index_to_grid(x, ind)
            #f <- function(i) series(x, vox[i,1], vox[i,2], vox[i,3])
            f <- function(i) series(x, ind[i])
            #lis <- lapply(seq_along(ind), function(i) f)
            #deferred_list(lis)
            deflist::deflist(f, length(ind))
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
            ##isplit <- split(1:length(clusters), clusters)


            #out <- vector(mode="list", length(isplit))

            #for (i in seq_along(isplit)) {
            #  out[[i]] <- function(i) {
            #    series_roi(x, isplit[[i]])
            #  }
            #}

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

#' @export
#' @rdname split_blocks-methods
setMethod(f="split_blocks", signature=signature(x="NeuroVec", indices="integer"),
          def = function(x, indices,...) {
            assert_that(length(indices) == dim(x)[4])
            isplit <- split(1:length(indices), indices)


            #out <- vector(mode="list", length(isplit))

            #for (i in seq_along(isplit)) {
            #  out[[i]] <- function(i) {
            #    sub_vector(x, isplit[[i]])
            #  }
            #}

            f <- function(i) {
              sub_vector(x, isplit[[i]])
            }

            #ret <- deferred_list(out)
            ret <- deflist::deflist(f, length(isplit), names=names(isplit))
            #names(ret) <- names(isplit)
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
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVol(dat, bspace)
          })



#' read_vec
#'
#' @description
#' This function loads an image volume from a file.
#'
#' @param file_name The name(s) of the file(s) to load. If more than one file_name is specified, 
#'   the files are loaded and concatenated.
#' @param indices The indices of the sub-volumes to load (e.g. if the file is 4-dimensional).
#' @param mask A mask defining the spatial elements to load.
#' @param mode The IO mode which is one of "normal", "mmap", "bigvec", or "filebacked".
#'
#' @return An \code{\linkS4class{NeuroVec}} object.
#'
#' @export
#' @note Memory-mapping a gzipped file is not currently allowed.
read_vec  <- function(file_name, indices=NULL, mask=NULL, mode=c("normal", "mmap", "bigvec", "filebacked")) {
  mode <- match.arg(mode)

  get_source <- function(x, indices, mask) {
    if (endsWith(x, ".nii")) {
      NeuroVecSource(x, indices, mask)
    } else if (endsWith(x, ".nii.gz")) {
      NeuroVecSource(x, indices, mask)
    } else if (endsWith(x, ".lv.h5")) {
      if (!is.null(indices)) {
        warning("do not support indices argument for .lv.h5 files.")
      }
      if (!is.null(mask)) {
        warning("do not support mask argument for .lv.h5 files.")
      }

      LatentNeuroVecSource(file_name)
    } else if (endsWith(x, ".h5")) {
      if (!is.null(indices)) {
        warning("do not support indices argument for .h5 files.")
      }
      if (!is.null(mask)) {
        warning("do not support mask argument for .h5 files.")
      }
      H5NeuroVecSource(file_name=x)
    }
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



#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="NeuroVec", y="NeuroVol"),
		def=function(x,y, ...) {
			.concat4D(x,y,...)
		})


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="NeuroVol", y="NeuroVec"),
  def=function(x,y, ...) {
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
#' @param drop whether to drop dimension of length 1
#' @rdname series-methods
setMethod(f="series", signature=signature(x="NeuroVec", i="integer"),
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



## @rdname series-methods
## @export
# setMethod("series", signature(x="NeuroVec", i="integer"),
#           def=function(x, i, j, k) {
#             if (missing(j) && missing(k)) {
#               vdim <- dim(x)[1:3]
#               mat <- arrayInd(i, vdim)
#             apply(mat, 1, function(i) x[i[1], i[2], i[3],])
#             } else {
#               x[i,j,k,]
#             }
#          })


#' @export
#' @param drop whether to drop dimension of length 1
#' @rdname series-methods
setMethod(f="series", signature=signature(x="H5NeuroVec", i="integer"),
          def=function(x,i,j,k, drop=TRUE) {
            if (missing(j) && missing(k)) {
              grid <- indexToGridCpp(i, dim(x)[1:3])
              callGeneric(x, grid)
            } else {
              ## could be solved with expand.grid, no?
              assert_that(length(i) == 1 && length(j) == 1 && length(k) ==1)
              ret <- x@obj[["data/elements"]][i,j,k,]
              if (drop) drop(ret) else ret
            }
          })


#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="numeric"),
		def=function(x, i, j, k) {
		  if (missing(j) && missing(k)) {
			  callGeneric(x,as.integer(i))
		  } else {
		    callGeneric(x,as.integer(i), as.integer(j), as.integer(k))
		  }
		})


#' @rdname series-methods
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
setMethod("series", signature(x="H5NeuroVec", i="numeric"),
          def=function(x, i, j, k) {
            if (missing(j) && missing(k)) {
              callGeneric(x,as.integer(i))
            } else {
              callGeneric(x,as.integer(i), as.integer(j), as.integer(k))
            }
          })



#' @rdname series-methods
#' @export
setMethod("series", signature(x="H5NeuroVec", i="matrix"),
          def=function(x,i) {
            assertthat::assert_that(ncol(i) == 3)

            d4 <- dim(x)[4]
            ir <- lapply(1:ncol(i), function(j) seq(min(i[,j]), max(i[,j])))

            ret <- x@obj[["data/elements"]][ir[[1]][1]:ir[[1]][[length(ir[[1]])]],
                                            ir[[2]][1]:ir[[2]][[length(ir[[2]])]],
                                            ir[[3]][1]:ir[[3]][[length(ir[[3]])]],,drop=FALSE]
            #ret2 <- apply(ret, 3L, c)
            ret2 <- t(array(ret, c(prod(dim(ret)[1:3]),dim(ret)[4])))
            if (ncol(ret2) != nrow(i)) {
              i2 <- apply(i, 2, function(ind) {
                ind - min(ind) + 1
              })

              i3 <- gridToIndex3DCpp(dim(ret)[1:3], i2)
              ret2[,i3,drop=FALSE]

            } else {
              ret2
            }
          })






#' @export
setAs(from="NeuroVec", to="matrix",
      function(from) {
        dm <- dim(from)
        d123 <- prod(dm[1:3])
        d4 <- dm[4]
        vals <- from[,]
        dim(vals) <- c(d123,d4)
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
#' @export
setMethod(f="as.matrix", signature=signature(x = "NeuroVec"), def=function(x) {
			as(x, "matrix")
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
#' @param nbit set nbit compression
#' @param compression compression level 1 to 9
#' @param chunk_dim the dimensions of each chunk
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


## NeuroVecSeq methods

#' Create a NeuroVecSeq Instance
#'
#' Constructs a NeuroVecSeq object to represent a variable-length list of NeuroVec objects.
#'
#' @param ... One or more instances of type \link{NeuroVec}.
#' @return A NeuroVecSeq object containing the provided NeuroVec objects, along with the associated space and length information.
#' @export
#'
#' @examples
#' v1 <- NeuroVec(array(0, c(5, 5, 5, 2)), space = NeuroSpace(dim = c(5, 5, 5, 2)))
#' v2 <- NeuroVec(array(0, c(5, 5, 5, 4)), space = NeuroSpace(dim = c(5, 5, 5, 4)))
#' v3 <- NeuroVec(array(0, c(5, 5, 5, 6)), space = NeuroSpace(dim = c(5, 5, 5, 6)))
#' vs <- NeuroVecSeq(v1, v2, v3)
#' blks <- split_blocks(vs, rep(1:3, each = 4))
NeuroVecSeq <- function(...) {
  vecs <- list(...)
  assert_that(all(map_lgl(vecs, ~ inherits(., "NeuroVec"))))
  sp <- space(vecs[[1]])
  lens <- map_dbl(vecs, function(x) dim(x)[4])
  sp <- add_dim(drop_dim(sp), sum(lens))

  new("NeuroVecSeq", space=sp, vecs=vecs, lens=lens)
}

#' @export
#' @rdname length-methods
setMethod("length", signature=c("NeuroVecSeq"),
          def=function(x) {
            sum(map_dbl(x@vecs, ~ length(.)))
          })


## @rdname series-methods
## @export
# setMethod("series", signature(x="NeuroVecSeq", i="integer"),
#           def=function(x, i, j, k) {
#             map(x@vecs, ~ series(., i,j,k)) %>% flatten_dbl()
#           })
#
#
## @rdname series-methods
## @export
# setMethod("series", signature(x="NeuroVecSeq", i="numeric"),
#           def=function(x, i, j, k) {
#             map(x@vecs, ~ series(., as.integer(i),as.integer(j),as.integer(k))) %>% flatten_dbl()
#           })
#
#
## @rdname series-methods
## @export
# setMethod("series", signature(x="NeuroVecSeq", i="matrix"),
#           def=function(x,i) {
#             do.call(rbind, map(x@vecs, ~ series(., i)))
#           })
#
#
## @rdname series-methods
## @export
# setMethod("series_roi", signature(x="NeuroVecSeq", i="matrix"),
#           def=function(x,i) {
#             rois <- map(x@vecs, ~ series_roi(., i))
#             if (length(rois) == 1) {
#               rois[[1]]
#             } else if (length(rois) == 2) {
#               concat(rois[[1]], rois[[2]])
#             } else {
#               f <- partial(concat, rois[[1]], rois[[2]])
#               do.call(f, rois[3:length(rois)])
#             }
#           })



#' [[
#' @rdname NeuroVecSeq-methods
#' @param x the object
#' @param i the indices
#' @export
setMethod(f="[[", signature=signature(x="NeuroVecSeq", i="numeric"),
          def = function(x, i) {

            assert_that(length(i) == 1 && i > 0 && i <= dim(x)[4])

            offsets <- cumsum(c(1, x@lens))[1:(length(x@lens))]
            vnum <- i - offsets
            vnum[vnum < 0] <- Inf
            bucket <- which.min(vnum)
            bucket_elnum <- vnum[bucket] + 1
            #print(bucket)
            #print(bucket_elnum)
            x@vecs[[bucket]][[bucket_elnum]]
          })

#' @noRd
setMethod(f="linear_access", signature=signature(x = "NeuroVecSeq", i = "numeric"),
          def = function (x, i) {
            ## inprog
            #browser()
            nels <- prod(dim(x)[1:3])
            els <- nels * x@lens
            csum <- cumsum(nels * x@lens) +1
            cels <- c(1, csum[-length(csum)])
            vnum <- find_seqnum(cels, i)
            offsets <- i - cels[vnum] + 1

            soff <- split(offsets, vnum)
            sind <- split(seq_along(vnum), vnum)

            res <- map(names(soff), function(vnum) {
              x@vecs[[as.integer(vnum)]][soff[[vnum]]]
            })

            out <- numeric(length(i))

            for (j in seq_along(res)) {
              out[sind[[j]]] <- res[[j]]
            }

            out

          })








