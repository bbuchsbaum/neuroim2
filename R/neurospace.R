#' @include all_class.R
{}
#' @include axis.R
{}


#' Constructor function for \code{\linkS4class{NeuroSpace}} class
#'
#' @param dim An integer vector describing the dimensions of the image grid.
#' @param origin A numeric vector representing the coordinate origin of the image space. If not provided, the default is set to a vector of zeroes with the same length as the dimensions.
#' @param spacing A numeric vector representing the real-valued voxel dimensions (e.g., in millimeters). If not provided, the default is set to a vector of ones with the same length as the dimensions.
#' @param axes An \code{\linkS4class{AxisSet}} object representing the image axes ordering. If not provided, the default axes are determined based on the NIFTI standard (Left-Posterior-Inferior).
#' @param trans A matrix representing the coordinate transformation associated with the image space. If not provided, the default is based on the NIFTI standard ("Neurological" orientation).
#'
#' @return An instance of the \code{\linkS4class{NeuroSpace}} class.
#'
#' @note Users rarely need to create a new \code{NeuroSpace} instance, as it will almost always be created automatically using information stored in an image header. If an existing image object is available, its \code{NeuroSpace} instance can be easily extracted with the \code{space} method.
#'
#' @examples
#' bspace <- NeuroSpace(c(64, 64, 64), origin = c(0, 0, 0), spacing = c(2, 2, 2))
#' bspace
#' origin(bspace)
#' axes(bspace)
#' trans(bspace)
#'
#' @export
#' @rdname neuro_space
NeuroSpace <- function(dim, spacing=NULL, origin=NULL, axes=NULL, trans=NULL) {

	if (is.null(spacing)) {
		spacing <- rep(1, min(length(dim), 3))
	}

	if (is.null(origin)) {
		origin <- rep(0, min(length(dim), 3))
	}

  if (length(origin) != length(spacing)) {
    stop("length of 'origin' must equal length of 'spacing'")
  }

  assertthat::assert_that(all(spacing > 0), msg="all dimensions must have sppacing > 0")
  assertthat::assert_that(all(dim > 0), msg="all dimensions must be > 0")

	if (is.null(trans)) {
		D <- min(length(dim), 3)
		trans <- diag(c(spacing,1))
		trans[1:D,D+1] <- origin
	}

  trans <- signif(trans,3)

	if (is.null(axes) && length(dim) >= 3) {
	  axes <- .nearestAnatomy(trans)
	} else if (is.null(axes) && length(dim) == 2) {
	  ### need .nearestAnatomy for 2d slice
	  ## TODO
	  axes <- AxisSet2D(LEFT_RIGHT, POST_ANT)
	}

	new("NeuroSpace", dim=as.integer(dim),
			origin=signif(origin,3),
			spacing=signif(spacing,3),
			axes=axes,
			trans=trans,
			inverse=signif(solve(trans),4))
}


#' show a \code{NeuroSpace}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("NeuroSpace"),
		def=function(object) {
			cat("NeuroSpace\n")
			cat("  Type           :", class(object), "\n")
			cat("  Dimension      :", object@dim, "\n")
			cat("  Spacing        :", paste(paste(object@spacing[1:(length(object@spacing)-1)], " X ", collapse=" "),
							object@spacing[length(object@spacing)], "\n"))
			cat("  Origin         :", paste(paste(object@origin[1:(length(object@origin)-1)], " X ", collapse=" "),
							object@origin[length(object@origin)], "\n"))
			#cat("  Axes           :", paste(object@axes@i@axis, object@axes@j@axis, object@axes@k@axis, "\n"))
			cat("  Coordinate Transform :", object@trans, "\n")


		}
)

#' add dimension to \code{\linkS4class{NeuroSpace}}
#' @export
#' @rdname add_dim-methods
setMethod(f="add_dim", signature=signature(x = "NeuroSpace", n="numeric"),
		def=function(x, n) {
			NeuroSpace(c(dim(x), n), origin=origin(x), spacing=spacing(x), axes=axes(x), trans=trans(x))
		})


#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x="NeuroSpace", dimnum="numeric"),
          def=function(x, dimnum) {
            D <- dim(x)
            stopifnot(length(D) >= 2)

            Dind <- seq(1,length(D))[-dimnum]
            if (ndim(x) > 3) {
              NeuroSpace(D[Dind], origin=origin(x)[Dind], spacing=spacing(x)[Dind],
                         axes=axes(x), trans=trans(x))
            } else {

              tx <- x@trans
              keep_col <- Dind
              keep_row <- which(apply(tx[,Dind], 1, function(x) !all(x==0)))
              tx <- rbind(cbind(tx[keep_row,keep_col], origin(x)[keep_row]), c(rep(0, length(Dind)), 1))

              NeuroSpace(D[-dimnum], origin=origin(x)[-dimnum], spacing=spacing(x)[-dimnum],
                         axes=drop_dim(axes(x), dimnum), trans=tx)
            }

          })

#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x = "NeuroSpace", dimnum="missing"),
		def=function(x) {

			D <- dim(x)
			stopifnot(length(D) >= 2)
			Dind <- 1:(length(D)-1)


			### TODO doesn't drop dimension in transformation matrix...
      ### brain vector's don't have th axis and these are incorrectly dropped
      if (ndim(x) > 3) {
			  NeuroSpace(D[Dind], origin=origin(x)[Dind], spacing=spacing(x)[Dind], axes=axes(x), trans=trans(x))
      } else {
        tx <- trans(x)
        tx <- rbind(cbind(tx[Dind,Dind], origin(x)[Dind]), c(rep(0, length(Dind)),1))
        NeuroSpace(D[Dind], origin=origin(x)[Dind], spacing=spacing(x)[Dind], axes=drop_dim(axes(x)), trans=tx)
      }
		})

#' dim
#'
#' @export
#' @param x the object
setMethod(f="dim", signature=signature(x = "NeuroSpace"),
		def=function(x) x@dim)


#' @export
#' @rdname ndim-methods
setMethod(f="ndim", signature=signature(x = "NeuroSpace"),
		def=function(x) length(x@dim))

#' @export
#' @rdname centroid-methods
setMethod(f="centroid", signature=signature(x = "NeuroSpace"),
          def=function(x)  {
            ind <- 1:prod(dim(x))
            colMeans(index_to_coord(x,ind))
          })


#' @export
#' @rdname dim_of-methods
setMethod(f="dim_of", signature=signature(x = "NeuroSpace", axis="NamedAxis"),
            function(x, axis) {
  dir <- abs(axis@direction)

  dnum <- which_dim(x,axis)
  dim(x)[dnum]
})

#' @export
#' @rdname which_dim-methods
setMethod(f="which_dim", signature=signature(x = "NeuroSpace", axis="NamedAxis"),
          function(x, axis) {
            dir <- abs(axis@direction)

            dnum = if(all(abs(x@axes@i@direction) == dir)) {
              1
            } else if (all(abs(x@axes@j@direction) == dir)) {
              2
            } else if (all(abs(x@axes@k@direction) == dir)) {
              3
            } else {
              stop(paste("cannot find matching axis of: ", axis))
            }

            dnum
          })


#' spacing
#'
#' @export
#' @rdname spacing-methods
setMethod(f="spacing", signature=signature(x = "NeuroSpace"), def=function(x) x@spacing)

#' bounds
#'
#' @export
#' @rdname bounds-methods
setMethod(f="bounds", signature=signature(x = "NeuroSpace"),
		def=function(x) {
		  c1 <- grid_to_coord(x, c(1,1,1))
		  c2 <- grid_to_coord(x, c(dim(x)[1], dim(x)[2], dim(x)[3]))
    	cbind(as.vector(c1),as.vector(c2))
		}
)


#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x="NeuroSpace", idx="numeric"),
          def=function(x, idx) {
            array.dim <- dim(x)
            .indexToGrid(idx, array.dim)
          })


#' @export
#' @rdname index_to_coord-methods
setMethod(f="index_to_coord", signature=signature(x="NeuroSpace", idx="numeric"),
          def=function(x, idx) {
            d <- min(3, ndim(x))
            grid <- index_to_grid(x, idx) - .5
            res <- trans(x) %*% t(cbind(grid[,1:d], rep(1,nrow(grid))))


            t(res[1:d,])
          })

#' @export
#' @rdname index_to_coord-methods
setMethod(f="index_to_coord", signature=signature(x="NeuroSpace", idx="integer"),
          def=function(x, idx) {
            grid <- index_to_grid(x, idx) - .5
            res <- trans(x) %*% t(cbind(grid, rep(1,nrow(grid))))
            t(res[1:ndim(x),])
          })

#' @export
#' @rdname index_to_coord-methods
setMethod(f="index_to_coord", signature=signature(x="NeuroVol", idx="integer"),
          def=function(x, idx) {
            callGeneric(space(x), as.numeric(idx))
          })


#' @export
#' @rdname index_to_coord-methods
setMethod(f="index_to_coord", signature=signature(x="NeuroVec", idx="integer"),
          def=function(x, idx) {
            callGeneric(space(x), as.numeric(idx))
          })


#' @export
#' @rdname coord_to_index-methods
setMethod(f="coord_to_index", signature=signature(x="NeuroSpace", coords="matrix"),
          def=function(x, coords) {
            grid = t(inverse_trans(x) %*% t(cbind(coords, rep(1, nrow(coords)))))
            grid_to_index(x, grid[,1:3] + .5)
          })

#' @export
#' @rdname coord_to_index-methods
setMethod(f="coord_to_index", signature=signature(x="NeuroSpace", coords="numeric"),
          def=function(x, coords) {
            coords <- matrix(coords, nrow=1)
            callGeneric(x,coords)
          })




#' @export
#' @rdname coord_to_grid-methods
setMethod(f="coord_to_grid", signature=signature(x="NeuroSpace", coords="matrix"),
          def=function(x, coords) {
            grid = t(inverse_trans(x) %*% t(cbind(coords, rep(1, nrow(coords)))))
            grid[,1:3] + 1
          })


#' @export
#' @rdname coord_to_grid-methods
setMethod(f="coord_to_grid", signature=signature(x="NeuroSpace", coords="numeric"),
          def=function(x, coords) {
            coords <- matrix(coords, nrow=1)
            callGeneric(x, coords)
          })


#' @export
#' @rdname grid_to_grid-methods
#' @importFrom purrr map_int
setMethod(f="grid_to_grid", signature=signature(x="NeuroSpace", vox="matrix"),
          def=function(x, vox) {

            nd <- ndim(x)
            stopifnot(ncol(vox) == nd)

            tx <- inverse_trans(x)[1:nd,1:nd]
            idx <- which(tx != 0, arr.ind=TRUE)
            tx[idx] <- 1 * sign(tx[idx])
            ovox <- tx %*% t(vox)
            offset <- map_dbl(1:nrow(tx), function(i) {
              if (any(tx[i,] < 0)) {
                dim(x)[i] + 1
              } else {
                0
              }

            })

            t(sweep(ovox, 1,offset, "+"))
          })

#' @export
#' @rdname grid_to_grid-methods
setMethod(f="grid_to_grid", signature=signature(x="matrix", vox="matrix"),
          def=function(x, vox) {
            nd <- ncol(x)-1
            stopifnot(ncol(vox) == nd)
            tx <- x[1:nd, 1:nd]
            ovox <- vox %*% tx[1:nd, 1:nd]
            offset <- map_int(1:ncol(tx), function(i) {
              if (any(tx[,i] < 0)) {
                dim(x)[i] + 1
              } else {
                0
              }
            })

            sweep(ovox, 2,offset, "+")

          })



#' @export
#' @rdname grid_to_coord-methods
setMethod(f="grid_to_coord", signature=signature(x="NeuroSpace", coords="matrix"),
          def=function(x, coords) {
            input <- t(cbind(coords-1, rep(1, nrow(coords))))
            ret <- t(trans(x) %*% input)
            md <- min(ndim(x), 3)
            ret[,1:md,drop=FALSE]
          })

#' @export
#' @rdname grid_to_coord-methods
setMethod(f="grid_to_coord", signature=signature(x="NeuroSpace", coords="numeric"),
          def=function(x, coords) {
            coords <- matrix(coords, ncol=length(coords))
            input <- t(cbind(coords-1, rep(1, nrow(coords))))
            ret <- t(trans(x) %*% input)

            md <- min(ndim(x), 3)
            ret[,1:md,drop=FALSE]

          })


#' @export
#' @rdname grid_to_coord-methods
setMethod(f="grid_to_coord", signature=signature(x="NeuroVol", coords="matrix"),
          def=function(x, coords) {
            callGeneric(space(x), coords)

          })



#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x="NeuroSpace", coords="matrix"),
		def=function(x, coords) {
			dx <- dim(x)

			if (ncol(coords) == 2) {
			  .gridToIndex(dim(x), coords)
			  #assert_that(ncol(coords) == 2)
			} else if (ncol(coords) == 3) {
			  assert_that(length(dx) >= 3)
			  .gridToIndex3D(dx[1:3], coords)
			} else if (ncol(coords) == 4 ){
			  .gridToIndex(dim(x), coords)
			} else {
			  stop(paste("grid_to_index: wrong dimensions: ndim = ", length(dim(x)), ", ncol(coords) = ", ncol(coords)))
			}
		})


#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x="NeuroSpace", coords="numeric"),
		def=function(x, coords) {
		  dx <- dim(x)
		  if (length(dx) == 2) {
		    assert_that(length(coords) == 2)
		    assert_that(coords[1] >= 1 && coords[1] <= dx[1] && coords[2] >= 1 && coords[2] <= dx[2])
		    (coords[2]-1)*dx[1] + coords[1]
		  } else {
		    assert_that(length(coords) == 3)
			  .gridToIndex3D(dim(x), matrix(coords, nrow=1, byrow=TRUE))
		  }
		}
)



#' @export
#' @rdname coord_to_index-methods
setMethod(f="coord_to_index", signature=signature(x="NeuroVol", coords="matrix"),
          def=function(x, coords) {
            assert_that(ncol(coords) == 3)
            callGeneric(space(x), coords)
          })



#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x="NeuroVec", idx="index"),
          def=function(x, idx) {
            callGeneric(space(x), idx)
          })

#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x="NeuroVec", idx="integer"),
          def=function(x, idx) {
            callGeneric(space(x), as.numeric(idx))
          })


#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x="NeuroVol", idx="index"),
          def=function(x, idx) {
            callGeneric(space(x), idx)
          })

#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid", signature=signature(x="NeuroVol", idx="integer"),
          def=function(x, idx) {
            callGeneric(space(x), as.numeric(idx))
          })

#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x="NeuroVol", coords="matrix"),
          def=function(x, coords) {
            assert_that(ncol(coords) == 3)
            array.dim <- dim(x)
            .gridToIndex3D(dim(x), coords)
          })


#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x="NeuroVol", coords="numeric"),
          def=function(x, coords) {
            assert_that(length(coords) == 3)
            array.dim <- dim(x)
            .gridToIndex3D(dim(x), matrix(coords, nrow=1, byrow=TRUE))
          })


#' @export
#' @rdname reorient-methods
setMethod(f="reorient", signature=signature(x = "NeuroSpace", orient="character"),
          def=function(x, orient) {

            anat <- findAnatomy3D(orient[1], orient[2], orient[3])
            pmat_orig <- perm_mat(x)
            pmat_new <- perm_mat(anat)


            tx <- t(pmat_new) %*% trans(x)[1:ndim(x),]
            tx <- rbind(tx,c(rep(0, ndim(x)),1))
            #itx <- zapsmall(MASS::ginv(tx))

            NeuroSpace(dim(x), spacing=spacing(x), axes=x@axes, trans=tx,
                          origin=tx[1:(ndim(x)) ,ndim(x)+1])

        }
)

#' @export
#' @rdname origin-methods
setMethod(f="origin", signature=signature(x = "NeuroSpace"), def=function(x) x@origin)

#' @export
#' @rdname origin-methods
setMethod(f="origin", signature=signature(x = "NeuroVol"), def=function(x) space(x)@origin)


#' @export
#' @rdname axes-methods
setMethod(f="axes", signature=signature(x = "NeuroSpace"), def=function(x) x@axes)


#' @export
#' @rdname trans-methods
setMethod(f="trans", signature=signature(x = "NeuroSpace"),def=function(x) x@trans)


#' @export
#' @rdname inverse_trans-methods
setMethod(f="inverse_trans", signature=signature(x = "NeuroSpace"), def=function(x) x@inverse)



#' @export
#' @rdname space-methods
setMethod(f="space", signature=signature(x = "NeuroSpace"), def=function(x) x)




