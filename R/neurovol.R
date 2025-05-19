#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseVector
#' @include all_class.R
#' @include all_generic.R
NULL


#' @importFrom methods new
#' @importFrom assertthat assert_that
NULL

#' NeuroVol: 3D Neuroimaging Volume Class
#'
#' @description
#' The \code{NeuroVol} class encapsulates 3D volumetric neuroimaging data. It provides methods
#' for accessing slices, performing spatial transformations, and integrating with the spatial
#' reference provided by \code{\linkS4class{NeuroSpace}}.
#'
#' @param data A 3D array containing the volumetric data.
#' @param space An object of class \code{\linkS4class{NeuroSpace}} defining the spatial properties.
#' @param label A character string providing a label for the volume (default: "").
#' @param indices An optional vector of indices for sparse representation (default: NULL).
#'
#' @return A \code{NeuroVol} object.
#'
#' @examples
#' bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
#' dat <- array(rnorm(64*64*64), c(64,64,64))
#' bvol <- NeuroVol(dat,bspace, label="test")
#' @export NeuroVol
#' @rdname NeuroVol
NeuroVol <- function(data, space,  label="", indices=NULL) {
	DenseNeuroVol(data,space, label=label, indices=indices)
}

#' Construct a DenseNeuroVol object
#'
#' Construct a \code{\linkS4class{DenseNeuroVol}} instance
#' @param data a three-dimensional \code{array}
#' @param space an instance of class \code{\linkS4class{NeuroSpace}}
#' @param label a \code{character} string
#' @param indices an optional 1-d index vector
#' @return \code{\linkS4class{DenseNeuroVol}} instance
#' @export DenseNeuroVol
#' @rdname DenseNeuroVol-class
DenseNeuroVol <- function(data, space, label="", indices=NULL) {

	if (length(dim(space)) != 3) {
		stop("DenseNeuroVol: space argument must have three dimensions")
	}

  if (is.matrix(data)) {
    if (nrow(data) == 1 || ncol(data) == 1) {
      data <- as.vector(data)
    }
  }

	if (length(data) == prod(dim(space)) && is.vector(data)) {
		dim(data) <- dim(space)
	}

	if (ndim(space) != 3) {
		stop("DenseNeuroVol: space argument must have three dimensions")
	}

	if (!all(dim(space) == dim(data))) {
		stop("DenseNeuroVol: data and space argument must have equal dimensions")
	}


	if (!is.null(indices)) {
	  assert_that(length(indices) == length(data))
		newdat <- array(0, dim(space))
		newdat[indices] <- data
		data <- newdat
	}


	new("DenseNeuroVol", data, space=space)

}






#' SparseNeuroVol
#'
#' Construct a \code{\linkS4class{SparseNeuroVol}} instance
#' @param data a numeric vector or ROIVol
#' @param space an instance of class \code{\linkS4class{NeuroSpace}}
#' @param indices a index vector indicating the 1-d coordinates of the data values
#' @param label a \code{character} string
#' @return \code{\linkS4class{SparseNeuroVol}} instance
#' @export SparseNeuroVol
#' @import Matrix
#' @details
#' Image data is backed by \code{Matrix::sparseVector}.
#' @examples
#' data <- 1:10
#' indices <- seq(1,1000, length.out=10)
#' bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
#' sparsevol <- SparseNeuroVol(data,bspace,indices=indices)
#' densevol <- NeuroVol(data,bspace,indices=indices)
#' sum(sparsevol) == sum(densevol)
#'
#' @rdname SparseNeuroVol-class
SparseNeuroVol <- function(data, space, indices=NULL, label="") {
  # Handle ROIVol input
  if (inherits(data, "ROIVol")) {
    if (is.null(indices)) {
      indices <- indices(data)
    }
    data <- as.numeric(data)
  }

  # Ensure data and indices are numeric vectors
  data <- as.numeric(data)
  indices <- as.numeric(indices)

  if (length(indices) != length(data)) {
    stop(paste("length of 'data' must equal length of 'indices'"))
  }

  # Create sparseVector with explicit class
  sv <- Matrix::sparseVector(x=data, i=indices, length=prod(dim(space)))

  new("SparseNeuroVol", data=sv, space=space)
}



#' LogicalNeuroVol
#'
#' This function constructs a \code{\linkS4class{LogicalNeuroVol}} instance.
#'
#' @param data A three-dimensional \code{array}, a 1D vector with length equal to \code{prod(dim(space))}, or a set of \code{indices} where elements are \code{TRUE}.
#' @param space An instance of class \code{\linkS4class{NeuroSpace}}.
#' @param label A \code{character} string.
#' @param indices An optional 1-d index vector.
#'
#' @return A \code{\linkS4class{LogicalNeuroVol}} instance.
#'
#' @examples
#' # Load an example brain mask
#' brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Convert the brain mask to a LogicalNeuroVol
#' logical_vol <- LogicalNeuroVol(brain_mask, space(brain_mask))
#'
#' @export LogicalNeuroVol
#' @rdname LogicalNeuroVol-class
LogicalNeuroVol <- function(data, space, label="", indices=NULL) {

	if (is.null(dim(data)) && length(data) == prod(dim(space))) {
		data <- array(data, dim(space))
	} else if (is.null(dim(data)) && !is.null(indices)) {
	  newdat <- array(FALSE, dim(space))
	  newdat[indices] <- data
	  data <- newdat
	}

	if (length(dim(data)) != 3) {
		stop("LogicalNeuroVol: data argument must have three dimensions")
	}

	if (ndim(space) != 3) {
		stop("LogicalVolume: space argument must have three dimensions")
	}


	if (!is.logical(data)) {
		D <- dim(data)
		data <- as.logical(data > 0)
		dim(data) <- D
	}

	new("LogicalNeuroVol", .Data=data, space=space)
}



#' Convert DenseNeuroVol to array
#'
#' This function converts a DenseNeuroVol object to an array.
#'
#' @param from A DenseNeuroVol object.
#'
#' @return An array resulting from the conversion.
#'
#' @keywords internal
#' @name as,DenseNeuroVol,array
setAs(from="DenseNeuroVol", to="array", def=function(from) from@.Data)






#' Convert SparseNeuroVol to array
#'
#' This function converts a SparseNeuroVol object to an array object.
#'
#' @param from A SparseNeuroVol object.
#'
#' @return An array object resulting from the conversion.
#'
#' @keywords internal
#' @name SparseNeuroVol,array
setAs(from="SparseNeuroVol", to="array", def=function(from) {
  vals <- as.numeric(from@data)
  array(vals, dim(from))
})


#' Convert SparseNeuroVol to numeric
#'
#' This function converts a SparseNeuroVol object to a numeric object.
#'
#' @param from A SparseNeuroVol object.
#'
#' @return A numeric object resulting from the conversion.
#' @keywords internal
#' @name SparseNeuroVol,numeric
setAs(from="SparseNeuroVol", to="numeric", def=function(from) {
  as.numeric(from@data)
})


#' Convert SparseNeuroVol to numeric
#'
#' @rdname as.numeric-methods
#' @param x the object to convert
#' @export
setMethod(f="as.numeric", signature=signature(x = "SparseNeuroVol"), def=function(x) {
  as(x, "numeric")
})



#' conversion from \code{\linkS4class{NeuroVol}} to \code{\linkS4class{LogicalNeuroVol}}
#'
#' @name as
#' @export
setAs(from="NeuroVol", to="LogicalNeuroVol", def=function(from) {
	LogicalNeuroVol(as.array(from), space(from))
})



#' conversion from \code{\linkS4class{NeuroVol}} to \code{array}
#'
#' @name as
#' @export
setAs(from="NeuroVol", to="array", def=function(from) from[,,])

#' Display NeuroVol Object
#'
#' @description Displays a formatted summary of a \code{NeuroVol} object.
#'
#' @param object A \code{NeuroVol} object.
#' @return Invisibly returns NULL after printing object information to the console.
#'
#' @importFrom crayon bold blue green red yellow silver
#' @importFrom utils object.size
#' @export
setMethod(f="show", signature=signature("NeuroVol"),
          def=function(object) {
            # Get space information and calculate stats
            sp <- space(object)
            val_range <- range(object, na.rm=TRUE)
            n_na <- sum(is.na(object))
            mem_size <- format(object.size(object), units="auto")
            total_voxels <- prod(dim(object))

            # Header
            cat("\n", crayon::bold(crayon::blue("=== NeuroVol Object ===")), "\n\n")

            # Basic Information
            cat(crayon::bold(crayon::yellow("* Basic Information")), "\n")
            cat("  ", crayon::silver("Type:"), "      ", class(object)[1], "\n", sep="")
            cat("  ", crayon::silver("Dimensions:"), " ",
                paste(dim(object), collapse=" x "),
                " (", crayon::green(mem_size), ")", "\n", sep="")
            cat("  ", crayon::silver("Total Voxels:"), " ",
                format(total_voxels, big.mark=","), "\n", sep="")

            # Data Properties
            cat("\n", crayon::bold(crayon::yellow("* Data Properties")), "\n", sep="")
            cat("  ", crayon::silver("Value Range:"), " [",
                crayon::blue(sprintf("%.2f", val_range[1])), ", ",
                crayon::blue(sprintf("%.2f", val_range[2])), "]", "\n", sep="")
            if (n_na > 0) {
                cat("  ", crayon::silver("Missing Values:"), " ",
                    crayon::red(format(n_na, big.mark=",")), " (",
                    sprintf("%.1f%%", 100*n_na/total_voxels), ")",
                    "\n", sep="")
            }

            # Spatial Properties
            cat("\n", crayon::bold(crayon::yellow("* Spatial Properties")), "\n", sep="")
            cat("  ", crayon::silver("Spacing:"), " ",
                paste(sprintf("%.2f", sp@spacing), collapse=" x "),
                crayon::silver(" mm"), "\n", sep="")
            cat("  ", crayon::silver("Origin:"), "  ",
                paste(sprintf("%.1f", sp@origin), collapse=", "),
                crayon::silver(" mm"), "\n", sep="")
            cat("  ", crayon::silver("Axes:"), "    ",
                crayon::green(sp@axes@i@axis), " x ",
                crayon::green(sp@axes@j@axis), " x ",
                crayon::green(sp@axes@k@axis), "\n", sep="")

            # Footer with usage hints
            cat(crayon::silver("\n======================================\n"))
            cat("\n", crayon::bold("Access Methods:"), "\n")
            cat(" ", crayon::silver("."), " Get Slice:  ",
                crayon::blue("slice(object, zlevel=10)"), "\n")
            cat(" ", crayon::silver("."), " Get Value:  ",
                crayon::blue("object[i, j, k]"), "\n")
            cat(" ", crayon::silver("."), " Plot:      ",
                crayon::blue("plot(object)"),
                crayon::silver(" # shows multiple slices"), "\n\n")
          })

#' show a \code{SparseNeuroVol}
#' @param object the object
#' @return Invisibly returns NULL after printing object information to the console.
#' @importFrom crayon bold blue green red yellow silver
#' @importFrom utils object.size
#' @importFrom Matrix which
#' @export
setMethod(f="show", signature=signature("SparseNeuroVol"),
          def=function(object) {
            # Get space information and calculate stats
            sp <- space(object)
            val_range <- range(as.numeric(object), na.rm=TRUE)
            n_na <- sum(is.na(as.numeric(object)))
            mem_size <- format(object.size(object), units="auto")
            total_voxels <- prod(dim(object))
            nonzero_voxels <- length(Matrix::which(object@data != 0))
            sparsity <- (total_voxels - nonzero_voxels) / total_voxels * 100

            # Header
            cat("\n", crayon::bold(crayon::blue("=== SparseNeuroVol Object ===")), "\n\n")

            # Basic Information
            cat(crayon::bold(crayon::yellow("* Basic Information")), "\n")
            cat("  ", crayon::silver("Type:"), "      ", class(object)[1], "\n", sep="")
            cat("  ", crayon::silver("Dimensions:"), " ",
                paste(dim(object), collapse=" x "),
                " (", crayon::green(mem_size), ")", "\n", sep="")

            # Sparsity Information
            cat("\n", crayon::bold(crayon::yellow("* Sparsity Properties")), "\n", sep="")
            cat("  ", crayon::silver("Total Voxels:  "),
                format(total_voxels, big.mark=","), "\n", sep="")
            cat("  ", crayon::silver("Active Voxels: "),
                crayon::green(format(nonzero_voxels, big.mark=",")),
                " (", sprintf("%.2f%%", 100-sparsity), ")", "\n", sep="")
            cat("  ", crayon::silver("Sparsity:     "),
                sprintf("%.2f%%", sparsity), "\n", sep="")

            # Data Properties
            cat("\n", crayon::bold(crayon::yellow("* Data Properties")), "\n", sep="")
            cat("  ", crayon::silver("Value Range:"), " [",
                crayon::blue(sprintf("%.2f", val_range[1])), ", ",
                crayon::blue(sprintf("%.2f", val_range[2])), "]", "\n", sep="")
            if (n_na > 0) {
                cat("  ", crayon::silver("Missing Values:"), " ",
                    crayon::red(format(n_na, big.mark=",")), " (",
                    sprintf("%.1f%%", 100*n_na/total_voxels), ")",
                    "\n", sep="")
            }

            # Spatial Properties
            cat("\n", crayon::bold(crayon::yellow("* Spatial Properties")), "\n", sep="")
            cat("  ", crayon::silver("Spacing:"), " ",
                paste(sprintf("%.2f", sp@spacing), collapse=" x "),
                crayon::silver(" mm"), "\n", sep="")
            cat("  ", crayon::silver("Origin:"), "  ",
                paste(sprintf("%.1f", sp@origin), collapse=", "),
                crayon::silver(" mm"), "\n", sep="")
            cat("  ", crayon::silver("Axes:"), "    ",
                crayon::green(sp@axes@i@axis), " x ",
                crayon::green(sp@axes@j@axis), " x ",
                crayon::green(sp@axes@k@axis), "\n", sep="")

            # Footer with usage hints
            cat(crayon::silver("\n======================================\n"))
            cat("\n", crayon::bold("Access Methods:"), "\n")
            cat(" ", crayon::silver("."), " Get Slice:     ",
                crayon::blue("slice(object, zlevel=10)"), "\n")
            cat(" ", crayon::silver("."), " Get Value:     ",
                crayon::blue("object[i, j, k]"), "\n")
            cat(" ", crayon::silver("."), " As Dense:      ",
                crayon::blue("as(object, \"DenseNeuroVol\")"), "\n")
            cat(" ", crayon::silver("."), " Active Indices: ",
                crayon::blue("which(object != 0)"), "\n\n")
          }
)



#' @rdname load_data-methods
#' @export
setMethod(f="load_data", signature=c(x="NeuroVolSource"),
		def=function(x) {

			meta <- x@meta_info
			nels <- prod(meta@dims[1:3])

			### for brain buckets, this offset needs to be precomputed ....
			offset <- (nels * (x@index-1)) * meta@bytes_per_element

			reader <- data_reader(meta, offset)
			dat <- read_elements(reader, nels)

			## bit of a hack to deal with scale factors
      if (.hasSlot(meta, "slope")) {

        if (meta@slope != 0) {
          dat <- dat*meta@slope
        }
      }

			close(reader)
			arr <- array(dat, meta@dims[1:3])

			bspace <- NeuroSpace(meta@dims[1:3], meta@spacing, meta@origin, meta@spatial_axes, trans(meta))
			DenseNeuroVol(arr, bspace, x)

		})

#' Constructor for NeuroVolSource
#'
#' @param input the input file name
#' @param index the image subvolume index
#' @export
#' @rdname NeuroVolSource-class
#' @return a new instance of type \code{NeuroVolSource}
NeuroVolSource <- function(input, index=1) {
	stopifnot(index >= 1)
	stopifnot(is.character(input))

	if (!file.exists(input)) {
		candidates <- Sys.glob(paste(input, "*", sep=""))
		if (length(candidates) > 0) {
			input <- candidates[1]
		}
	}

	stopifnot(file.exists(input))

	meta_info <- read_header(input)

	if (length(meta_info@dims) < 4) {
		if (index > 1) {
			stop(sprintf("index cannot be greater than 1 for a 3D image (dimensions: %s)",
			            paste(meta_info@dims, collapse="x")))
		}
	} else if (index > meta_info@dims[4]) {
		stop(sprintf("index %d exceeds available volumes in 4D image (max: %d)",
		            index, meta_info@dims[4]))
	}

	new("NeuroVolSource", meta_info=meta_info, index=as.integer(index))
}

#' Load an image volume from a file
#'
#' @param file_name the name of the file to load
#' @param index the index of the volume (e.g. if the file is 4-dimensional)
#' @return an instance of the class \code{\linkS4class{DenseNeuroVol}}
#'
#' @examples
#' fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
#' x <- read_vol(fname)
#' print(dim(x))
#' space(x)
#'
#' @export read_vol
read_vol  <- function(file_name, index=1) {
	src <- NeuroVolSource(file_name, index)
	load_data(src)
}


#' @export
#' @rdname slices-methods
#' @param x A NeuroVol object
#' @return A deflist object containing functions that return 2D slices of the volume along the z-axis.
#'         The length of the deflist equals the number of slices in the z dimension.
setMethod(f="slices", signature=signature(x="NeuroVol"),
          def = function(x) {
            nslices <- dim(x)[3]
            f <- function(i) slice(x, i, 3)
            #lis <- lapply(1:nslices, function(i) f)
            #deferred_list(lis)
            deflist::deflist(f, nslices)
          })


#' @export
setMethod(f="[", signature=signature(x = "NeuroVol", i = "ROICoords", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            callGeneric(x, i@coords)
          }
)



#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="DenseNeuroVol", y="missing"),
          def=function(x,y,...) {
            x
          })


#' @note dimensions of x and y must be equal
#' @export
#' @rdname concat-methods
setMethod(f="concat", signature=signature(x="DenseNeuroVol", y="DenseNeuroVol"),
		def=function(x,y,...) {
			.concat4D(x,y,...)
		})


#' @export
#' @rdname map_values-methods
setMethod(f="map_values", signature=signature(x="NeuroVol", lookup="list"),
          def=function(x,lookup) {
            if (is.null(names(lookup))) {
              names(lookup) <- seq_along(lookup)
            }
            out <- match(x,as.numeric(names(lookup)))
            DenseNeuroVol(unlist(lookup[out]), space(x))
          })

#' @export
#' @rdname map_values-methods
setMethod(f="map_values", signature=signature(x="NeuroVol", lookup="matrix"),
          def=function(x,lookup) {
            if (ncol(lookup) != 2) {
              stop("map_values: lookup matrix must have two columns: column 1 is key, column 2 is value")
            } else if (nrow(lookup) < 1) {
              stop("map_values: lookup matrix have at least one row")
            }

            m <- match(as.vector(x), lookup[,1])
            outv <- lookup[m,2]
            outv[is.na(outv)] <- 0
            NeuroVol(array(outv, dim(x)), space(x))
          })


#' @export
#' @rdname split_fill-methods
setMethod(f="split_fill", signature=signature(x="NeuroVol", fac="factor", FUN="function"),
		def=function(x,fac,FUN) {
			stopifnot(length(x) == length(fac))
			S <- split(1:length(x), fac, drop=TRUE)
			res <- map(S, function(ind) {
						X <- FUN(x[ind])
						if (length(X) == length(ind)) {
							X
						} else {
							rep(X, length(ind))
						}
			})

			ovol <- x
			ovol[1:length(x)] <- unsplit(res, fac)
			ovol

		})

#' @keywords internal
#' @noRd
fast.expand.grid <- function(seq1,seq2, constant) {
  cbind(Var1 = rep.int(seq1, length(seq2)),
        Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))),
        Var3=constant)
}

#' @export
#' @rdname slice-methods
setMethod(f="slice", signature=signature(x="NeuroVol", zlevel="numeric", along="numeric",
                                         orientation="missing"),
          def=function(x, zlevel, along, orientation) {
            stopifnot(along >= 1 && along <=3)
            stopifnot(zlevel >= 1 && zlevel <= dim(x)[along])

            imslice <- switch(as.character(along),
                              "1"=x[zlevel,,],
                              "2"=x[,zlevel,],
                              "3"=x[,,zlevel])

            NeuroSlice(imslice, drop_dim(space(x), along))


          })

#' @export
#' @rdname slice-methods
setMethod(f="slice", signature=signature(x="NeuroVol", zlevel="numeric", along="NeuroSpace",
                                         orientation="AxisSet3D"),
          def=function(x, zlevel, along, orientation) {


            xdim <- dim_of(along, orientation@i)
            ydim <- dim_of(along, orientation@j)
            zdim <- dim_of(along, orientation@k)


            stopifnot(zlevel >= 1 && zlevel <= zdim)

            vox <- as.matrix(fast.expand.grid(seq(1,xdim), seq(1,ydim), zlevel))

            gg <- grid_to_grid(along, vox)

            imslice <- x[gg]

            NeuroSlice(matrix(imslice, xdim,ydim), drop_dim(along,which_dim(along, orientation@k)))


          })


#' @title Convert NeuroVol to a mask
#' @description This method converts a NeuroVol object to a mask by setting all positive values to TRUE and all non-positive values to FALSE.
#' @param x A NeuroVol object to convert to a mask.
#' @param indices A missing argument; not used in this method.
#' @return A LogicalNeuroVol object representing the mask created from the input NeuroVol.
#' @export
#' @rdname as.mask-methods
setMethod(f="as.mask", signature=signature(x="NeuroVol", indices="missing"),
          def=function(x) {
            LogicalNeuroVol(x > 0, space(x))
          })

#' @title Convert NeuroVol to a mask with specified indices
#' @description This method converts a NeuroVol object to a mask by setting the specified indices to TRUE and the remaining elements to FALSE.
#' @param x A NeuroVol object to convert to a mask.
#' @param indices A numeric vector containing the indices of the input NeuroVol that should be set to TRUE in the resulting mask.
#' @return A LogicalNeuroVol object representing the mask created from the input NeuroVol with specified indices.
#' @export
#' @rdname as.mask-methods
setMethod(f="as.mask", signature=signature(x="NeuroVol", indices="numeric"),
          def=function(x, indices) {
            M <- array(0, dim(x))
            M[indices] <- 1
            LogicalNeuroVol(M, space(x))
          })


#' @export
#' @rdname coord_to_grid-methods
setMethod(f="coord_to_grid", signature=signature(x="NeuroVol", coords="matrix"),
          def=function(x, coords) {
            callGeneric(space(x), coords)
          })

#' @export
#' @rdname coord_to_grid-methods
setMethod(f="coord_to_grid", signature=signature(x="NeuroVol", coords="numeric"),
          def=function(x, coords) {
            coords <- matrix(coords, nrow=1)
            callGeneric(x, coords)
          })


#' @importFrom dbscan kNN
#' @keywords internal
#' @noRd
.pruneCoords <- function(coord.set,  vals,  mindist=10) {

	if (NROW(coord.set) == 1) {
		1
	}

	.prune <- function(keepIndices) {
		if (length(keepIndices) <= 2) {
			keepIndices
		} else {
		  ret <- try(dbscan::kNN(coord.set[keepIndices,], coord.set[keepIndices,], k=2))

			#ret <- rflann::Neighbour(coord.set[keepIndices,], coord.set[keepIndices,], k=2, build="kdtree", cores=0, checks=1)
			#ind <- ret$indices[, 2]
			#ds <- sqrt(ret$distances[, 2])

		  ind <- ret$id[,2]
		  ds <- ret$distances[,2]

			v <- vals[keepIndices]
			ovals <- v[ind]
			pruneSet <- ifelse(ds < mindist & ovals > v,  TRUE, FALSE)

			if (any(pruneSet)) {
				Recall(keepIndices[!pruneSet])
			} else {
				keepIndices
			}
		}


	}

	.prune(1:NROW(coord.set))



}


#' @title Create a patch set from a NeuroVol object
#' @description This function creates a patch set from a NeuroVol object given specified dimensions
#' @return A deferred list of patches.
#' @param x a NeuroVol object
#' @param dims the dimensions of the patch
#' @param ... additional args
#' @export
#' @rdname patch_set-methods
setMethod(f="patch_set", signature=signature(x="NeuroVol",
                                             dims="numeric",
                                             mask="missing"),
          def=function(x, dims, ...) {
            mask <- LogicalNeuroVol(array(1, dim(x)), space=space(x))
            callGeneric(x, dims, mask)
          })


#' @title Create a patch set from a NeuroVol object
#' @description This function creates a patch set from a NeuroVol object given specified dimensions and a mask.
#' @return A deferred list of patches.
#' @param x a NeuroVol object
#' @param dims the dimensions of the patch
#' @param mask the mask defining the valid patch centers
#' @param ... additional args
#' @export
#' @rdname patch_set-methods
setMethod(f="patch_set", signature=signature(x="NeuroVol",
                                            dims="numeric",
                                            mask="LogicalNeuroVol"),
          def=function(x, dims, mask, ...) {
            if (!all(dims > 0)) {
              stop("all 'dims' must be greater than zero")
            }

            if (!all(dims %% 2 == 1)) {
              stop("all 'dims' must be odd numbers")
            }

            template <- as.matrix(expand.grid
                                    (x=seq(seq(ceiling(-dims[1]/2),floor(dims[1]/2))),
                                    y=seq(seq(ceiling(-dims[2]/2),floor(dims[2]/2))),
                                    z=seq(seq(ceiling(-dims[3]/2),floor(dims[3]/2)))))

            xdim <- dim(x)
            grid <- index_to_grid(mask, which(mask>0))

            f <- function(i) {
              g <- grid[i,]
              m <- t(t(template) + g)
              m[,1] <- pmax(pmin(m[,1], xdim[1]), 1)
              m[,2] <- pmax(pmin(m[,2], xdim[2]), 1)
              m[,3] <- pmax(pmin(m[,3], xdim[3]), 1)
              ret <- x[m]
              attr(ret, "idx") <- grid_to_index(x, m)
              ret

            }

            patches <- deflist::deflist(f, nrow(grid))
            patches
            #patches <- deferred_list(lapply(1:nrow(grid), function(i) f))

          })


#' @rdname values-methods
#' @export
setMethod(f="values", signature=signature(x="DenseNeuroVol"),
          def=function(x) {
            x@.Data
          })

#' @rdname values-methods
#' @export
setMethod(f="values", signature=signature(x="SparseNeuroVol"),
          def=function(x) {
            as.vector(x@data)
          })

#' apply a kernel function to a \code{\linkS4class{NeuroVol}}
#'
#' @rdname map-methods
#' @param mask restrict application of kernel to masked area
#' @export
setMethod(f="mapf", signature=signature(x="NeuroVol", m="Kernel"),
          def=function(x, m, mask=NULL) {
            ovol <- array(0, dim(x))
            hwidth <- map_dbl(m@width, function(d) ceiling(d/2 -1)) + 1
            xdim <- dim(x)[1]
            ydim <- dim(x)[2]
            zdim <- dim(x)[3]

            if (!is.null(mask)) {
              if (!all.equal(dim(mask), dim(ovol))) {
                stop(paste("mask must have same dimensions as input volume"))
              }

              grid <- index_to_grid(mask, which(mask != 0))
            } else {
              grid <- as.matrix(expand.grid(i=hwidth[1]:(xdim - hwidth[1]),
                                            j=hwidth[2]:(ydim - hwidth[2]),
                                            k=hwidth[3]:(zdim - hwidth[3])))
            }

            res <- apply(grid, 1, function(vox) {
              loc <- t(t(m@voxels) +vox)
              ivals <- x[loc]
              sum(ivals * m@weights)
            })

            ovol[grid] <- res
            NeuroVol(ovol, space(x))
          })


#' @export
#' @importFrom utils capture.output
#' @rdname show-methods
setMethod("show", "Kernel", function(object) {
  # Get dimensions of kernel
  kernel_dims <- dim(object@weights)

  # Create header string
  header <- sprintf("Kernel object of dimensions: %s", paste(kernel_dims, collapse=" x "))

  # Get kernel width info
  width_str <- sprintf("Kernel widths: %s", paste(object@width, collapse=" x "))

  # Format weights matrix/array nicely
  weights_str <- capture.output(print(round(object@weights, 4)))
  weights_str <- paste(weights_str, collapse="\n")

  # Print everything
  cat(header, "\n")
  cat(width_str, "\n")
  cat("\nKernel weights:\n")
  cat(weights_str, "\n")
})




#' Find connected components in NeuroVol
#'
#' @export
#' @importFrom purrr map_int
#' @param threshold threshold defining lower intensity bound for image mask
#' @param cluster_table return cluster_table
#' @param local_maxima return table of local maxima
#' @param local_maxima_dist the distance used to define minum distance between local maxima
#' @rdname conn_comp-methods
setMethod(f="conn_comp", signature=signature(x="NeuroVol"),
	def=function(x, threshold=0, cluster_table=TRUE, local_maxima=TRUE, local_maxima_dist=15,...) {
		mask <- (x > threshold)
		stopifnot(any(mask))

		comps <- conn_comp_3D(mask,...)

		grid <- as.data.frame(index_to_grid(mask, which(mask>0)))
		colnames(grid) <- c("x", "y", "z")
		locations <- split(grid, comps$index[comps$index>0])

		ret <- list(index=ClusteredNeuroVol(mask, clusters=comps$index[mask>0]), size=NeuroVol(comps$size, space(x)), voxels=locations)

		if (cluster_table) {
			maxima <- do.call(rbind, lapply(locations, function(loc) {
				if (nrow(loc) == 1) {
					loc
				} else {
					vals <- x[as.matrix(loc)]
					loc[which.max(vals),]
				}
			}))
			N <- comps$size[as.matrix(maxima)]
			Area <- N * prod(spacing(x))
			maxvals <- x[as.matrix(maxima)]
			ret$cluster_table <- data.frame(index=1:NROW(maxima),
			                                x=maxima[,1], y=maxima[,2], z=maxima[,3],
			                                N=N, Area=Area, value=maxvals)
		}

		if (local_maxima) {
			#if (all(map_int(locations, NROW) == 1)) {
      #
			#}
			coord.sets <- lapply(locations, function(loc) {
				sweep(as.matrix(loc), 2, spacing(x), "*")
			})

			loc.max <- do.call(rbind, mapply(function(cset, i) {
				idx <- .pruneCoords(as.matrix(cset), x[as.matrix(locations[[i]])], mindist=local_maxima_dist)
				maxvox <- as.matrix(locations[[i]])[idx,,drop=F]
				cbind(i, maxvox)
			}, coord.sets, 1:length(coord.sets), SIMPLIFY=FALSE))


			loc.max <- cbind(loc.max, x[loc.max[, 2:4, drop=F]])

			row.names(loc.max) <- 1:NROW(loc.max)
			colnames(loc.max) <- c("index", "x", "y", "z", "value")
			ret$local_maxima <- loc.max
		}

		ret


	})




### TODO when the source voulme is an AFNI BRIK the output file does not preserve orientation
### e.g. write_vol(x, space(BRIK), ...)
### if BRIK was RAI, output NIFTI file reverts to LPI


#' @export
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="NeuroVol", file_name="character", format="missing", data_type="missing"),
		def=function(x, file_name) {
			write_nifti_volume(x, file_name)
		})


#' @export
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="ClusteredNeuroVol", file_name="character", format="missing", data_type="missing"),
          def=function(x, file_name) {
            callGeneric(as(x, "DenseNeuroVol"), file_name)
          })




#' @export
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="NeuroVol", file_name="character", format="character", data_type="missing"),
		def=function(x, file_name, format) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, file_name)
			} else {
				stop(paste("sorry, cannot write format: ", format))
			}
		})


#' @export
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="ROIVol", file_name="character", format="character", data_type="missing"),
          def=function(x, file_name, format) {
            if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
              callGeneric(as.dense(x), file_name)
            } else {
              stop(paste("sorry, cannot write format: ", format))
            }
          })



#' @export write_vol
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="NeuroVol", file_name="character", format="missing", data_type="character"),
		def=function(x, file_name, data_type) {
			write_nifti_volume(x, file_name, data_type)

		})

#' as.logical
#'
#' Convert NeuroVol to \code{\linkS4class{LogicalNeuroVol}}
#'
#' the image values will be converted to using R base function \code{as.logical} and wrapped in \code{LogicalNeuroVol}
#'
#' @param x the object
#' @return an instance of \code{\linkS4class{LogicalNeuroVol}}
#' @rdname as.logical-methods
#' @export
setMethod(f="as.logical", signature=signature(x = "NeuroVol"), def=function(x) {
			vals <- as.logical(as.vector(x))
			LogicalNeuroVol(vals, space(x))
})

#' @rdname as.sparse-methods
#' @export
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVol", mask="LogicalNeuroVol"),
          def=function(x, mask) {
            assert_that(all(dim(x) == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))
            dat <- x[mask]
            bvec <- SparseNeuroVol(data=dat, space=space(x),indices=which(mask>0))

})


#' @rdname partition-methods
setMethod("partition", signature=signature(x="LogicalNeuroVol", k="integer"),
          def=function(x,k) {
            assert_that(k>1)
            idx <- which(x != 0)
            cds <- index_to_coord(x,idx)
            kres <- kmeans(cds, centers=k, iter.max=200)
            ClusteredNeuroVol(x, kres$cluster)
          })

#' @rdname partition-methods
setMethod("partition", signature=signature(x="LogicalNeuroVol", k="numeric"),
          def=function(x,k) {
            callGeneric(x, as.integer(k))
          })

#' @rdname partition-methods
setMethod("partition", signature=signature(x="DenseNeuroVol", k="numeric"),
          def=function(x,k) {
            callGeneric(as.logical(x), as.integer(k))
          })

#' @rdname as.sparse-methods
#' @export
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVol", mask="numeric"),
          def=function(x, mask) {
            m <- as.integer(mask)
            bvol <- SparseNeuroVol(x[m], space(x), indices=m)
            bvol
          })




#' @rdname linear_access-methods
#' @export
setMethod(f="linear_access", signature=signature(x = "SparseNeuroVol", i = "numeric"),
          def=function (x, i) {
            x@data[as.numeric(i)]
          })




#' @export
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "numeric", j = "numeric", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            if (missing(k)) {
              k <- seq(1, dim(x)[3])
            }

            nel <- length(i) * length(j) * length(k)
            mind <- cbind(rep(i, length.out=nel), rep(j, each=length(i)), rep(k, each=length(i) * length(j)))

            ind <- grid_to_index(x, mind)
            dat <- x@data[as.numeric(ind)]
            array(dat, c(length(i), length(j), length(k)))

          })





#' plot an NeuroVol as a series of 2D slices
#'
#' @rdname plot-methods
#' @param x the object to display
#' @param cmap a color map consisting of a vector of colors in hex format (e.g. \code{gray(n=255)})
#' @param zlevels the series of slice indices to display.
#' @param along the axis along which to slice
#' @param irange the intensity range indicating the low and high values of the color scale.
#' @param thresh a 2-element vector indicating the lower and upper transparency thresholds.
#' @param alpha the level of alpha transparency
#' @param bgvol a background volume that serves as an image underlay (currently ignored).
#' @param bgcmap a color map for backround layer consisting of a vector of colors in hex format (e.g. \code{gray(n=255)})
#' @export
#' @importFrom graphics plot
#' @import colorplane
#' @examples
#'
#' dat <- matrix(rnorm(100*100), 100, 100)
#' slice <- NeuroSlice(dat, NeuroSpace(c(100,100)))
#' \donttest{
#' plot(slice)
#' }
#'
#' vol_data <- array(rnorm(10 * 10 * 10), c(10, 10, 10))
#' vol <- NeuroVol(vol_data, NeuroSpace(c(10, 10, 10)))
#' \donttest{
#' plot(vol, along=3)
#' }
setMethod("plot", signature=signature(x="NeuroVol"),
          def=function(x,
                       cmap=gray(seq(0,1,length.out=255)),
                       zlevels=unique(round(seq(1, dim(x)[3], length.out=6))),
                       along=3,
                       irange=range(x, na.rm=TRUE),
                       thresh=c(0,0),
                       alpha=1,
                       bgvol=NULL,
                       bgcmap=gray(seq(0,1,length.out=255))) {

            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            if (!is.null(bgvol)) {
              assert_that(all(dim(x) == dim(bgvol)))
              assert_that(all(spacing(x) == spacing(bgvol)))
            }

            # Create a data frame of all the slices specified in zlevels
            df1 <- do.call(rbind, purrr::map(zlevels, function(i) {
              if (!is.null(bgvol)) {
                bgslice <- slice(bgvol, zlevel=i, along=along)
                bgplane <- colorplane::IntensityColorPlane(as.numeric(bgslice), cols=bgcmap)
                bgcols <- colorplane::map_colors(bgplane)
              }

              imslice <- slice(x, zlevel=i, along=along)
              implane <- colorplane::IntensityColorPlane(as.numeric(imslice), cols=cmap, alpha=alpha)
              fgcols <- colorplane::map_colors(implane, threshold=thresh, irange=irange)

              if (!is.null(bgvol)) {
                fgcols <- colorplane::as_hexcol(colorplane::blend_colors(bgcols, fgcols, alpha=alpha))
              } else {
                fgcols <- colorplane::as_hexcol(fgcols)
              }

              cds <- index_to_coord(space(imslice), 1:length(imslice))
              data.frame(x=cds[,1], y=cds[,2], z=i, value=as.vector(fgcols))
            }))

            {y = value = NULL} # to appease R CMD check

            p <- ggplot2::ggplot(df1, ggplot2::aes(x=x, y=y)) +
              ggplot2::coord_fixed() +
              ggplot2::geom_raster(ggplot2::aes(fill=value)) +
              ggplot2::scale_fill_identity() +
              ggplot2::facet_wrap(~ z, labeller = ggplot2::labeller(z = function(z) paste("Slice:", z))) +
              ggplot2::ggtitle("Brain Slices") +
              ggplot2::theme_void() +
              ggplot2::theme(
                strip.background = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(face="bold", size=10),
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "cm")
              )

            p
          })


