#' @include all_class.R
NULL
#' @include all_generic.R
NULL



#' NeuroVec
#'
#' constructor function for virtual class \code{\linkS4class{NeuroVec}}
#'
#' @param data the image data which can be a \code{matrix}, a 4d \code{array}, or a list of \code{NeuroVols}.
#'        If the latter, the geometric space of the data \code{NeuroSpace} will be inferred from the constituent volumes,
#'        which must all be identical.
#' @param space a \code{\linkS4class{NeuroSpace}} object. Does not ned to be included if \code{data} argument is a list of \code{NeuroVols}
#' @param mask an optional \code{array} of type \code{logical}
#' @param label a label of type \code{character}
#' @return a concrete instance of \code{\linkS4class{NeuroVec}} class.
#' If \code{mask} is provided then \code{\linkS4class{SparseNeuroVec}}, otherwise \code{\linkS4class{DenseNeuroVec}}
#' @export NeuroVec
#' @rdname NeuroVec-class
NeuroVec <- function(data, space=NULL, mask=NULL, label="") {
  if (is.list(data)) {
    space <- space(data[[1]])
    space <- add_dim(space, length(data))
    data <- do.call(cbind, lapply(data, function(x) as.vector(x)))

  }

	if (prod(dim(space)) != length(data)) {
		stop("dimensions of data argument do not match dimensions of space argument")
	}
	if (is.null(mask)) {
		DenseNeuroVec(data,space, label)
	} else {
		SparseNeuroVec(data,space,mask,label)
	}

}


#' DenseNeuroVec
#'
#' constructor function for class \code{\linkS4class{DenseNeuroVec}}
#'
#' @param data a 4-dimensional \code{array} or a 2-dimension \code{matrix} that is either nvoxels by ntime-points or ntime-points by nvoxels
#' @param space a \code{\linkS4class{NeuroSpace}} object
#' @param label a label of type \code{character}
#' @return \code{\linkS4class{DenseNeuroVec}} instance
#' @export DenseNeuroVec
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
	}



	new("DenseNeuroVec", .Data=data,  space=space)

}




#' load_data
#' @return an instance of class \code{\linkS4class{NeuroVec}}
#' @importFrom RNifti readNifti
#' @rdname load_data-methods
setMethod(f="load_data", signature=c("NeuroVecSource"),
		def=function(x) {

			meta <- x@meta_info

			stopifnot(length(meta@dims) == 4)

			nels <- prod(meta@dims[1:4])
			ind <- x@indices

			## use RNifti, fails to work with other formats, though...
			arr <- RNifti::readNifti(meta@data_file)

			## bit of a hack to deal with scale factors
			if (.hasSlot(meta, "slope")) {

        if (meta@slope != 0) {
			    arr <- arr* meta@slope
        }
			}

      bspace <- NeuroSpace(c(meta@dims[1:3], length(ind)),meta@spacing, meta@origin,
                           meta@spatial_axes, trans(meta))
			DenseNeuroVec(arr[,,,ind,drop=FALSE], bspace, x)

		})



#' NeuroVecSource
#'
#' Construct a \code{\linkS4class{NeuroVecSource}} object
#'
#' @param file_name name of the 4-dimensional image file
#' @param indices the subset of integer volume indices to load -- if \code{NULL} then all volumes will be loaded
#' @param mask image volume indicating the subset of voxels that will be loaded. If provided, function returns \code{\linkS4class{SparseNeuroVecSource}}
#' @return a instance deriving from \code{\linkS4class{NeuroVecSource}}
#'
#' @details If a \code{mask} is supplied then it should be a \code{\linkS4class{LogicalNeuroVol}} or \code{\linkS4class{NeuroVol}} instance. If the latter, then the mask will be defined by nonzero elements of the volume.
#'
#' @rdname NeuroVecSource
#' @importFrom assertthat assert_that
#' @export
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
		SparseNeuroVecSource(meta_info, as.integer(indices), mask)
	}

}


#' Get length of \code{NeuroVec}. This is the numbe rof volumes in the volume vector (e.g. the 4th image dimension)
#'
#' @export
#' @rdname length-methods
setMethod("length", signature=c("NeuroVec"),
		def=function(x) {
			dim(x)[4]
		})



#' read_vol_list
#'
#' load a list of image volumes and return a \code{\linkS4class{NeuroVec}} instance
#'
#' @param file_names a list of files to load
#' @param mask an optional mask indicating subset of voxels to load
#' @return an instance of class \code{\linkS4class{NeuroVec}}
#' @export read_vol_list
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


		SparseNeuroVec(mat[mask,], dspace, mask=mask, label=map_chr(meta_info, function(m) m@label))

	}
}


setAs("DenseNeuroVec", "array", function(from) from@.Data)

setAs("NeuroVec", "array", function(from) from[,,,])

#' show a \code{NeuroVecSource}
#' @param object the object
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




#' show a \code{NeuroVec}
#' @param object the object
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
setMethod(f="sub_vector", signature=signature(x="DenseNeuroVec", i="numeric"),
          def=function(x, i) {
            assertthat::assert_that(max(i) <= dim(x)[4])
            xs <- space(x)
            dat <- x[,,,i]

            newdim <- c(dim(x)[1:3], length(i))
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVec(dat, bspace)
          })


#' @export
#' @rdname vols-methods
setMethod(f="vols", signature=signature(x="NeuroVec", indices="numeric"),
          def = function(x, indices) {
            assert_that(min(indices) > 0 && max(indices) <= dim(x)[4])
            lis <- lapply(indices, function(i) function(i) x[[i]])
            deferred_list(lis)
          })

#' @export
#' @rdname vols-methods
setMethod(f="vols", signature=signature(x="NeuroVec", indices="missing"),
          def = function(x) {
            f <- function(i) x[[i]]
            lis <- lapply(1:(dim(x)[4]), function(i) f)
            deferred_list(lis)
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="missing"),
          def = function(x) {
            ind <- 1:prod(dim(x)[1:3])
            vox <- index_to_grid(x, ind)
            f <- function(i) x[vox[i,1], vox[i,2], vox[i,3],]
            lis <- lapply(ind, function(i) f)
            deferred_list(lis)
          })

#' @export
#' @rdname vectors-methods
setMethod(f="vectors", signature=signature(x="NeuroVec", subset="numeric"),
          def = function(x,subset) {
            ind <- subset
            assert_that(max(ind) < prod(dim(x)[1:3]))
            vox <- index_to_grid(x, ind)
            f <- function(i) x[vox[i,1], vox[i,2], vox[i,3],]
            lis <- lapply(seq_along(ind), function(i) f)
            deferred_list(lis)
          })


#' [[
#' @rdname NeuroVec-methods
#' @param i the volume index
#' @export
setMethod(f="[[", signature=signature(x="NeuroVec", i="numeric"),
          def = function(x, i) {
            xs <- space(x)
            dat <- x[,,,i]
            newdim <- dim(x)[1:3]
            bspace <- NeuroSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            DenseNeuroVol(dat, bspace)
          })



#' read_vec
#'
#' load an image volume from a file
#'
#' @param file_name the name of the file to load
#' @param indices the indices of the sub-volumes to load (e.g. if the file is 4-dimensional)
#' @param mask a mask defining the spatial elements to load
#' @return an \code{\linkS4class{NeuroVec}} object
#' @export
read_vec  <- function(file_name, indices=NULL, mask=NULL) {
	src <- NeuroVecSource(file_name, indices, mask)
	load_data(src)
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
            grid <- coords(i)
            callGeneric(x, grid)
          })


#' @rdname series-methods
#' @export
setMethod("series_roi", signature(x="NeuroVec", i="ROICoords"),
          def=function(x,i) {
            rvol <- series(x, i)
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
setMethod("series_roi", signature(x="NeuroVec", i="LogicalNeuroVol"),
          def=function(x,i) {
            mat <- as.matrix(series(x, i))
            ROIVec(space(x), coords=index_to_grid(i, which(i == TRUE)), data=as.matrix(mat))
          })

#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroVec", i="numeric"),
		def=function(x, i, j, k) {
			if (missing(j) && missing(k)) {
				vdim <- dim(x)[1:3]
				mat <- arrayInd(i, vdim)
				apply(mat, 1, function(i) x[i[1], i[2], i[3],])
			} else {
				x[i,j,k,]
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
              stop("series_roi: must provide either 1D 'i' or 3D ('i', 'k', 'j') vector indices")
            }
            else {
              vox <- cbind(i,j,k)
              callGeneric(x, as.matrix(vox))
            }


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
  out = list()
  for (i in 1:dim(x)[4]) {
    out[[i]] <- takeVolume(x,i)
  }

  out
})


#' convert a \code{DenseNeuroVec} to a matrix
#'
#' @rdname as.matrix-methods
#' @param x the object
#' @export
setMethod(f="as.matrix", signature=signature(x = "DenseNeuroVec"), def=function(x) {
			as(x, "matrix")
		})


#' @rdname as.sparse-methods
#' @export
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVec", mask="LogicalNeuroVol"),
          def=function(x, mask) {
            assert_that(all(dim(x)[1:3] == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))

            vdim <- dim(x)[1:3]
            dat <- as.matrix(x)[mask == TRUE,]
            bvec <- SparseNeuroVec(dat, space(x), mask)

          })


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
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="missing"),
		def=function(x, file_name) {
			write_nifti_vector(x, file_name)
		})


#' @export
#' @rdname write_vec-methods
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="character", data_type="missing"),
		def=function(x, file_name, format) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, file_name)
			} else {
				stop(paste("sorry, cannot write format: ", format))
			}
		})


#' @export write_vec
#' @rdname write_vec-methods
#' @aliases write_vec,NeuroVec,character,missing,character,ANY-method
setMethod(f="write_vec",signature=signature(x="NeuroVec", file_name="character", format="missing", data_type="character"),
		def=function(x, file_name, data_type) {
			write_nifti_vector(x, file_name, data_type)

		})


