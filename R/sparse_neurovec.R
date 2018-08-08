#' @include all_class.R
{}
#' @include all_generic.R
{}



#' SparseNeuroVecSource
#'
#' constructs a SparseNeuroVecSource object
#'
#' @param meta_info an object of class \code{\linkS4class{MetaInfo}}
#' @param indices a vector of 1D indices
#' @param mask a 3D \code{array} of type \code{logical}
#' @export
#' @rdname SparseNeuroVecSource-class
SparseNeuroVecSource <- function(meta_info, indices, mask) {


	stopifnot(length(dim(meta_info)) >= 3)
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


#' SparseNeuroVec
#'
#' constructs a SparseNeuroVec object
#'
#' @param data an array which can be a \code{matrix} or 4-D \code{array}
#' @param space a NeuroSpace instance
#' @param mask a 3D \code{array} of type \code{logical}
#' @param source the data source -- an instance of class \code{\linkS4class{FileSource}}
#' @param label associated sub-image labels
#' @export
#' @examples
#'
#' bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
#' mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
#' mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
#' svec <- SparseNeuroVec(mat, bspace,mask)
#' length(indices(svec)) == sum(mask)
#' @rdname SparseNeuroVec-class
SparseNeuroVec <- function(data, space, mask, source=NULL, label="") {
	stopifnot(inherits(space, "NeuroSpace"))

	if (!inherits(mask, "LogicalNeuroVol")) {
		mspace <- NeuroSpace(dim(space)[1:3], spacing(space), origin(space), axes(space), trans(space))
		mask <- LogicalNeuroVol(as.logical(mask), mspace)
	}

	stopifnot(inherits(mask, "LogicalNeuroVol"))


	D4 <- if (is.matrix(data)) {
		Nind <- sum(mask == TRUE)
		if (nrow(data) == Nind) {
			data <- t(data)
			nrow(data)
		} else if (ncol(data) == Nind) {
			nrow(data)
		} else {
			stop(paste("matrix with dim:", dim(data), " does not match mask cardinality: ", Nind))
		}
	} else if (length(dim(data)) == 4) {
		mat <- apply(data, 4, function(vals) vals)
		data <- t(mat[mask==TRUE,])
		dim(data)[4]
	}

	if (ndim(space) == 3) {
		space <- add_dim(space, nrow(data))
	}

  stopifnot(ndim(space) == 4)

	new("SparseNeuroVec", space=space, mask=mask,
	    data=data, map=IndexLookupVol(space(mask), as.integer(which(mask))))

}




#' @export
#' @rdname load_data-methods
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

#' @export
#' @rdname indices-methods
setMethod(f="indices", signature=signature(x="SparseNeuroVec"),
          def=function(x) {
            indices(x@map)
          })


#' @export
#' @rdname coords-methods
setMethod(f="coords", signature=signature(x="SparseNeuroVec"),
          def=function(x,i) {
            if (missing(i)) {
              return(coords(x@map, indices(x@map)))
            }
            coords(x@map, i)
          })


#' @export
#' @rdname series-methods
setMethod(f="series", signature=signature(x="SparseNeuroVec", i="matrix"),
         def=function(x,i) {
           idx <- grid_to_index(x@mask, i)
           callGeneric(x,idx)
         })



setMethod("series", signature(x="SparseNeuroVec", i="numeric"),
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
 setMethod("series", signature(x="SparseNeuroVec", i="integer"),
		 def=function(x,i, j, k) {
			 if (missing(j) && missing(k)) {
				 idx <- lookup(x, as.integer(i))
				 idx.nz <- idx[idx!=0]
				 if (length(idx.nz) == 0) {
					 matrix(0, dim(x)[4], length(i))
				 } else {
					 mat <- matrix(0, dim(x)[4], length(i))
					 mat[, idx !=0] <- x@data[,idx.nz]
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

#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="SparseNeuroVec", y="missing"),
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


#' @export
#' @rdname lookup-methods
setMethod(f="lookup", signature=signature(x="SparseNeuroVec", i="numeric"),
         def=function(x,i) {
            lookup(x@map, i)
          })


setMethod(f="linear_access", signature=signature(x = "SparseNeuroVec", i = "numeric"),
          def=function (x, i) {
            nels <- prod(dim(x)[1:3])
            n <- as.integer(i / nels) + 1
            offset <- i %% nels

            ll <- lookup(x, offset)
            nz <- which(ll > 0)
            idx2d <- cbind(n[nz], ll[nz])
            vals <- x@data[idx2d]

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
setMethod(f="[", signature=signature(x = "SparseNeuroVec", i = "numeric", j = "numeric"),
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
            oval[rep(keep, length(m))] <- x@data[indmat]

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




#' @export
setAs(from="SparseNeuroVec", to="matrix",
		  function(from) {
		    ind <- indices(from)
		    out <- matrix(0, dim(from)[4], length(ind))
		    out[, ind] <- from@data
		    out
			  #from@data
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
            cat("an instance of class",  class(object), "\n\n")
            cat("   dimensions: ", dim(object), "\n")
            cat("   voxel spacing: ", spacing(object), "\n")
            cat("   cardinality: ", length(object@map@indices))
            cat("\n\n")

          })




