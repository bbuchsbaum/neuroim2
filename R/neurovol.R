#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseVector
#' @importFrom rflann Neighbour
#' @include all_class.R
#' @include all_generic.R
NULL

#' make_vol
#'
#' Construct a \code{\linkS4class{NeuroVol}} instance, using default (dense) implementation
#'
#' @param data an optional one- or three-dimensional \code{vector} or \code{array}
#' @param refvol an instance of class \code{\linkS4class{NeuroVol}} containing the reference space for the new volume.
#' @param label an optional \code{character} string
#' @param indices an optional 1d vector of indices in to the 3d space
#' @return \code{\linkS4class{DenseNeuroVol}} instance
#' @examples
#' bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
#' dat <- array(rnorm(64*64*64), c(64,64,64))
#' bvol <- NeuroVol(dat,bspace, label="test")
#' bvol2 <- make_vol(dat, bvol)
#' all.equal(as.array(bvol),as.array(bvol2))
#' data <- 1:10
#' indices = seq(1,1000, length.out=10)
#' bvol3 <- make_vol(data,bvol,indices=indices)
#' sum(bvol3) == sum(data)
#' @export make_vol
make_vol <- function(data=NULL, refvol, label="", indices=NULL) {
  if (is.null(data)) {
	  DenseNeuroVol(array(0, dim(refvol)),space(refvol),label,indices)
  } else {
    DenseNeuroVol(data,space(refvol),label,indices)
  }
}

#' NeuroVol
#'
#' Construct a \code{\linkS4class{NeuroVol}} instance, using default (dense) implementation
#' @param data a three-dimensional \code{array}
#' @param space an instance of class \code{\linkS4class{NeuroSpace}}
#' @param label a \code{character} string to identify volume
#' @param indices an 1D vector that gives the linear indices of the associated \code{data} vector
#' @return a \code{\linkS4class{DenseNeuroVol}} instance
#' @examples
#' bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
#' dat <- array(rnorm(64*64*64), c(64,64,64))
#' bvol <- NeuroVol(dat,bspace, label="test")
#' print(bvol)
#' @export NeuroVol
#' @rdname NeuroVol
NeuroVol <- function(data, space,  label="", indices=NULL) {
	DenseNeuroVol(data,space, label=label, indices=indices)
}

#' DenseNeuroVol
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
		newdat <- array(0, dim(space))
		newdat[indices] <- data
		data <- newdat
	}


	new("DenseNeuroVol", data, space=space)

}






#' SparseNeuroVol
#'
#' Construct a \code{\linkS4class{SparseNeuroVol}} instance
#' @param data a numeric vector
#' @param space an instance of class \code{\linkS4class{NeuroSpace}}
#' @param indices a index vector indicating the 1-d coordinates of the data values
#' @param label a \code{character} string
#' @return \code{\linkS4class{SparseNeuroVol}} instance
#' @export SparseNeuroVol
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
#'
#'
#' @rdname SparseNeuroVol-class
SparseNeuroVol <- function(data, space, indices=NULL, label="") {
  if (length(indices) != length(data)) {
    stop(paste("length of 'data' must equal length of 'indices'"))
  }

  sv <- Matrix::sparseVector(x=data, i=indices, length=prod(dim(space)))

  new("SparseNeuroVol", data=sv, space=space)
}



#' LogicalNeuroVol
#'
#' Construct a \code{\linkS4class{LogicalNeuroVol}} instance
#' @param data a three-dimensional \code{array}, a 1D vector with length equal to \code{prod(dim(space))}, or a set of \code{indices} where elements are \code{TRUE}
#' @param space an instance of class \code{\linkS4class{NeuroSpace}}
#' @param label a \code{character} string
#' @param indices an optional 1-d index vector
#' @return \code{\linkS4class{LogicalNeuroVol}} instance
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
		data <- as.logical(data)
		dim(data) <- D
	}

	new("LogicalNeuroVol", .Data=data, space=space)
}



#' conversion from DenseNeuroVol to array
#' @rdname as-methods
#' @name as
setAs(from="DenseNeuroVol", to="array", def=function(from) from@.Data)



#' conversion from \code{SparseNeuroVol} to \code{array}
#' @rdname as-methods
#' @name as
setAs(from="SparseNeuroVol", to="array", def=function(from) {
  vals <- as.numeric(from@data)
  array(vals, dim(from))
})


#' conversion from SparseNeuroVol to numeric
#' @rdname as-methods
#' @name as
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


#' conversion from NeuroVol to LogicalNeuroVol
#' @rdname as-methods
#' @name as
setAs(from="NeuroVol", to="LogicalNeuroVol", def=function(from) {
	LogicalNeuroVol(as.array(from), space(from))
})

#' conversion from DenseNeuroVol to LogicalNeuroVol
#' @name as
#' @rdname as-methods
setAs(from="DenseNeuroVol", to="LogicalNeuroVol", def=function(from) {
	LogicalNeuroVol(as.array(from), space(from))
})



#' conversion from NeuroVol to array
#' @rdname as-methods
#' @name as
setAs(from="NeuroVol", to="array", def=function(from) from[,,])

#' show a \code{NeuroVol}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("NeuroVol"),
          def=function(object) {
            sp <- space(object)
            cat("NeuroVol\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(signif(sp@spacing[1:(length(sp@spacing)-1)],2), " X ", collapse=" "),
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(signif(sp@origin[1:(length(sp@origin)-1)],2), " X ", collapse=" "),
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")

          }
)

#' show a \code{SparseNeuroVol}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("SparseNeuroVol"),
          def=function(object) {
            sp <- space(object)
            cat("NeuroVol\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(signif(sp@spacing[1:(length(sp@spacing)-1)],2), " X ", collapse=" "),
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(signif(sp@origin[1:(length(sp@origin)-1)],2), " X ", collapse=" "),
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")
            cat("  Cardinality    :", length(which(object@data>0)))
          }
)



#' load a NeuroVol
#' @export load_data
#' @rdname load_data-methods
## TODO reduce code duplication with load_data#NeuroVecSource
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

	if ( length(meta_info@dims) < 4 && index > 1) {
		stop("index cannot be greater than 1 for a image with dimensions < 4")
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
#' fname <- system.file("extdata", "global_mask.nii", package="neuroim")
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
setMethod(f="slices", signature=signature(x="NeuroVol"),
          def = function(x) {
            nslices <- dim(x)[3]
            f <- function(i) slice(x, i, 3)
            lis <- lapply(1:nslices, function(i) f)
            deferred_list(lis)
          })

setMethod(f="[", signature=signature(x = "NeuroVol", i = "ROICoords", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            browser()
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
#' @export concat
#' @rdname concat-methods
setMethod(f="concat", signature=signature(x="DenseNeuroVol", y="DenseNeuroVol"),
		def=function(x,y,...) {
			.concat4D(x,y,...)
		})


#' @export
#' @rdname map_values-methods
setMethod(f="map_values", signature=signature(x="NeuroVol", lookup="list"),
          def=function(x,lookup) {
            out <- lookup[unlist(x)]
            DenseNeuroVol(unlist(out), space(x))
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


#' @export split_fill
#' @rdname split_fill-methods
setMethod(f="split_fill", signature=signature(x="NeuroVol", fac="factor", FUN="function"),
		def=function(x,fac,FUN) {
			stopifnot(length(x) == length(fac))
			S <- split(1:length(x), fac, drop=TRUE)
			res <- sapply(S, function(ind) {
						X <- FUN(x[ind])
						if (length(X) == length(ind)) {
							X
						} else {
							rep(X, length(ind))
						}
					}, simplify=FALSE)

			ovol <- x
			ovol[1:length(x)] <- unsplit(res, fac)
			ovol

		})

#' @keywords internal
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

            #browser()

            stopifnot(zlevel >= 1 && zlevel <= zdim)

            vox <- as.matrix(fast.expand.grid(seq(1,xdim), seq(1,ydim), zlevel))

            gg <- gridToGrid(along, vox)

            imslice <- x[gg]

            NeuroSlice(matrix(imslice, xdim,ydim), drop_dim(along,which_dim(along, orientation@k)))


          })







#' @export
#' @rdname as.mask-methods
setMethod(f="as.mask", signature=signature(x="NeuroVol", indices="missing"),
          def=function(x) {
            LogicalNeuroVol(x > 0, space(x))
          })

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


#' @importFrom rflann Neighbour
#' @keywords internal
.pruneCoords <- function(coord.set,  vals,  mindist=10) {

	if (NROW(coord.set) == 1) {
		1
	}

	.prune <- function(keepIndices) {
		if (length(keepIndices) == 1) {
			keepIndices
		} else {
			ret <- rflann::Neighbour(coord.set[keepIndices,], coord.set[keepIndices,], k=2)
			ind <- ret$indices[, 2]
			ds <- sqrt(ret$distances[, 2])
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

#' @rdname patch_set-methods
#' @export
setMethod(f="patch_set", signature=signature(x="NeuroVol",
                                             dims="numeric",
                                             mask="missing"),
          def=function(x, dims, ...) {
            mask <- LogicalNeuroVol(array(1, dim(x)), space=space(x))
            callGeneric(x, dims, mask)

          })


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

            template <- as.matrix(expand.grid(x=seq(seq(ceiling(-dims[1]/2),floor(dims[1]/2))),
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
              x[m]

            }

            patches <- deferred_list(lapply(1:nrow(grid), function(i) f))

          })




#' apply a kernel function to a \code{\linkS4class{NeuroVol}}
#'
#' @rdname map-methods
#' @param mask restrict application of kernel to masked area
#' @export
setMethod(f="map", signature=signature(x="NeuroVol", m="Kernel"),
          def=function(x, m, mask=NULL) {
            ovol <- array(0, dim(x))
            hwidth <- sapply(m@width, function(d) ceiling(d/2 -1)) + 1
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
              loc <- sweep(m@voxels, 2, vox, "+")
              ivals <- x[loc]
              if (all(ivals == 0)) {
                0
              } else {
                sum(ivals * m@weights)
              }
            })

            ovol[grid] <- res
            NeuroVol(ovol, space(x))
          })







#' find connected components in NeuroVol
#'
#' @export
#' @param threshold threshold defining lower intensity bound for image mask
#' @param cluster_table return cluster_table
#' @param local_maxima return table of local maxima
#' @param local_maxima_dist the distance used to define minum distance between local maxima
#' @rdname conn_comp-methods
setMethod(f="conn_comp", signature=signature(x="NeuroVol"),
	def=function(x, threshold=0, cluster_table=TRUE, local_maxima=TRUE, local_maxima_dist=15) {
		mask <- (x > threshold)
		stopifnot(any(mask))

		comps <- connComp3D(mask)

		grid <- as.data.frame(index_to_grid(mask, which(mask>0)))
		colnames(grid) <- c("x", "y", "z")
		locations <- split(grid, comps$index[comps$index>0])

		ret <- list(size=NeuroVol(comps$size, space(x)), index=NeuroVol(comps$index, space(x)), voxels=locations)



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
			ret$cluster_table <- data.frame(index=1:NROW(maxima), x=maxima[,1], y=maxima[,2], z=maxima[,3], N=N, Area=Area, value=maxvals)
		}

		if (local_maxima) {
			if (all(sapply(locations, NROW) == 1)) {

			}
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


#' @export write_vol
#' @rdname write_vol-methods
setMethod(f="write_vol",signature=signature(x="NeuroVol", file_name="character", format="missing", data_type="character"),
		def=function(x, file_name, data_type) {
			write_nifti_volume(x, file_name, data_type)

		})

#' as.logical
#'
#' Convert NeuroVol to \code{linkS4class{LogicalNeuroVol}}
#'
#' the image values will be converted to using R base function \code{as.logical} and wrapped in \code{LogicalNeuroVol}
#'
#' @param x the object
#' @return an instance of \code{linkS4class{LogicalNeuroVol}}
#' @rdname as.logical-methods
#' @export
setMethod(f="as.logical", signature=signature(x = "NeuroVol"), def=function(x) {
			vals <- as.logical(as.vector(x))
			LogicalNeuroVol(vals, space(x))
})

#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVol", mask="LogicalNeuroVol"),
          def=function(x, mask) {
            assert_that(all(dim(x) == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))

            dat <- x[mask]
            bvec <- SparseNeuroVec(dat, space(x))

})

#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseNeuroVol", mask="numeric"),
          def=function(x, mask) {
            m <- as.integer(mask)
            bvec <- SparseNeuroVol(x[m], space(x), indices=m)

          })



#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "numeric", j = "numeric", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            if (missing(k)) {
              k <- seq(1, dim(x)[3])
            }

            n <- length(i) * length(j) * length(k)
            mind <- cbind(rep(i, length.out=n), rep(j, each=length(i)), rep(k, each=length(i) * length(j)))

            ind <- gridToIndex(x, mind)
            x@data[ind]

        })


#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop drop dimension
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "numeric", j = "missing", drop="missing"),
          def=function (x, i, j, k, ..., drop) {
            if (missing(k) && nargs() == 4) {
              x@data[i]
            } else {
              callGeneric(x, i=i,  seq(1,dim(x)[2]), k, drop)
            }
         }
)

#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            ind <- gridToIndex(x,i)
            x@data[ind]
          }
)

#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "missing", j = "missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            if (missing(k)) {
              x@data
            } else {
              callGeneric(x, seq(1, dim(x)[1]), seq(1, dim(x)[2]), k, ...)
            }
          }
)

#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseNeuroVol", i = "missing", j = "numeric", drop="ANY"),
          def=function (x, i, j, k,  ..., drop=TRUE) {
            callGeneric(x, seq(1, dim(x)[1]), j, k,...)
          }
)




