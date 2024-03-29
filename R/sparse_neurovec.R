#' @include all_class.R
{}
#' @include all_generic.R
{}



#' SparseNeuroVecSource
#'
#' constructs a SparseNeuroVecSource object
#'
#' @param meta_info an object of class \code{\linkS4class{MetaInfo}}
#' @param indices an optional vector of 1D indices the subset of volumes to load
#' @param mask a logical 3D \code{array},  a logical 1D \code{vector} or a \code{LogicalNeuroVol}
#' @rdname SparseNeuroVecSource-class
#' @keywords internal
#' @noRd
SparseNeuroVecSource <- function(meta_info, indices=NULL, mask) {

  if (is.null(indices)) {
    indices <- seq(1, dim(meta_info)[4])
  }

	assert_that(length(dim(meta_info)) >= 3)
	stopifnot(all(indices >= 1 & indices <= dim(meta_info)[4]))

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
		stop("illegal mask argument with dim: ", paste(dim(mask), collapse=", "))
	}

  if (!inherits(mask, "LogicalNeuroVol")) {
    mspace <- NeuroSpace(dim(mask),  meta_info@spacing, meta_info@origin, meta_info@spatial_axes)
    mask <- LogicalNeuroVol(mask, mspace)
  }

	stopifnot(all(dim(mask) == D))

	new("SparseNeuroVecSource", meta_info=meta_info, indices=indices, mask=mask)
}


#' @keywords internal
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

  stopifnot(inherits(mask, "LogicalNeuroVol"))

  D4 <- if (is.matrix(data)) {
    Nind <- sum(mask == TRUE)
    if (nrow(data) == Nind) {
      data <- t(data)
      assert_that(ncol(data) == cardinality, msg="data matrix must must match cardinality of `mask`")
      nrow(data)
    } else if (ncol(data) == Nind) {
      assert_that(ncol(data) == cardinality, msg="data matrix must must match cardinality of `mask`")
      nrow(data)
    } else {
      stop(paste(
        "matrix with dim:",
        dim(data),
        " does not match mask cardinality: ",
        Nind
      ))
    }
  } else if (length(dim(data)) == 4) {
    mat <- apply(data, 4, function(vals)
      vals)
    data <- t(mat[mask == TRUE, ])
    dim(data)[4]
  }

  if (ndim(space) == 3) {
    space <- add_dim(space, D4)
  }

  stopifnot(ndim(space) == 4)

  list(mask=mask, data=data, space=space)
}


#' Construct a SparseNeuroVec Object
#'
#' Constructs a SparseNeuroVec object for efficient representation and manipulation
#' of sparse neuroimaging data with many zero or missing values.
#'
#' @param data A matrix or a 4-D array containing the neuroimaging data. The dimensions of the data should be consistent with the dimensions of the provided NeuroSpace object and mask.
#' @param space A \link{NeuroSpace} object representing the dimensions and voxel spacing of the neuroimaging data.
#' @param mask A 3D array, 1D vector of type logical, or an instance of type \link{LogicalNeuroVol}, which specifies the locations of the non-zero values in the data.
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
SparseNeuroVec <- function(data, space, mask) {
	stopifnot(inherits(space, "NeuroSpace"))
  p <- prep_sparsenvec(data,space, mask)

	new("SparseNeuroVec", space=p$space, mask=p$mask,
	    map=IndexLookupVol(space(p$mask), as.integer(which(p$mask))), data=p$data)

}





#' @keywords internal
#' @noRd
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


 #' @export
 #' @rdname series-methods
 #' @param j index for 2nd dimension
 #' @param k index for 3rd dimension
 setMethod("series", signature(x="AbstractSparseNeuroVec", i="integer"),
		 def=function(x,i, j, k) {
			 if (missing(j) && missing(k)) {
				 idx <- lookup(x, as.integer(i))
				 idx.nz <- idx[idx!=0]
				 if (length(idx.nz) == 0) {
					 matrix(0, dim(x)[4], length(i))
				 } else {
					 mat <- matrix(0, dim(x)[4], length(i))
					 #mat[, idx !=0] <- x@data[,idx.nz]
					 #browser()
					 mat[, idx !=0] <- matricized_access(x, idx.nz)
					 mat
				 }
			 } else {
				 vdim <- dim(x)
				 idx <- gridToIndex3DCpp(vdim[1:3], cbind(i,j,k))
				 #slicedim <- vdim[1] * vdim[2]
				 #idx <- slicedim*(k-1) + (j-1)*vdim[1] + i
				 callGeneric(x, idx)
			 }

		 })


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

#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "matrix"),
          def=function (x, i) {
            x@data[i]
          })

#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "integer"),
          def=function (x, i) {
            x@data[,i]
          })

#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "SparseNeuroVec", i = "numeric"),
          def=function (x, i) {
            x@data[,i]
          })


#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "matrix"),
          def=function (x, i) {
            x@data[i]
          })

#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "integer"),
          def=function (x, i) {
            x@data[,i]
          })

#' @noRd
#' @keywords internal
setMethod(f="matricized_access", signature=signature(x = "BigNeuroVec", i = "numeric"),
          def=function (x, i) {
            x@data[,i]
          })

#' @noRd
#' @keywords internal
setMethod(f="linear_access", signature=signature(x = "AbstractSparseNeuroVec", i = "numeric"),
          def=function (x, i) {
            nels <- prod(dim(x)[1:3])
            n <- ceiling(i/nels)
            offset <- i %% nels
            offset[offset == 0] <- nels

            ll <- lookup(x, offset)
            nz <- which(ll > 0)

            #if (length(nz) == 0) {
            #  return(numeric(length(i)))
            #}

            idx2d <- cbind(n[nz], ll[nz])
            vals <- matricized_access(x, idx2d)
            ##vals <- x@data[idx2d]

            ovals <- numeric(length(i))
            ovals[nz] <- vals
            ovals
          })



#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
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

#' @export
#' @rdname sub_vector-methods
setMethod(f="sub_vector", signature=signature(x="SparseNeuroVec", i="numeric"),
          def=function(x, i) {
            idx <- which(x@mask > 0)
            bspace <- drop_dim(space(x))

            res <- lapply(i, function(i) x@data[i,])
            res <- do.call("cbind", res)
            SparseNeuroVec(res, bspace, x@mask)
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


#' as.matrix
#'
#' convert SparseNeuroVec to matrix
#' @rdname as.matrix-methods
#' @export
setMethod(f="as.matrix", signature=signature(x = "SparseNeuroVec"), def=function(x) {
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
#' @export
setMethod("show",
          signature=signature(object="SparseNeuroVec"),
          def=function(object) {
            cat(class(object), "\n\n")
            cat("   Dimension: ", dim(object), "\n")
            cat("   Spacing: ", spacing(object), "\n")
            cat("   Origin: ", origin(space(object)), "\n")
            cat("   Cardinality: ", length(object@map@indices))
            cat("\n\n")
          })




