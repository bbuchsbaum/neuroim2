#' @include all_class.R
{}
#' @include axis.R
{}

#' Spatial Reference System for Neuroimaging Data
#'
#' @title NeuroSpace: Spatial Reference System for Neuroimaging Data
#'
#' @description
#' The \code{NeuroSpace} class defines the spatial properties and coordinate system of
#' neuroimaging data. It encapsulates all information needed to map between voxel indices
#' and real-world coordinates, including dimensions, voxel spacing, origin, axis orientation,
#' and coordinate transformations.
#'
#' @param dim An integer vector specifying the dimensions of the image grid. Must be positive.
#' @param spacing A numeric vector specifying the physical size of each voxel (typically in
#'   millimeters). Must be positive. If NULL, defaults to ones.
#' @param origin A numeric vector specifying the real-world coordinates of the first voxel.
#'   If NULL, defaults to zeros.
#' @param axes An \code{\linkS4class{AxisSet}} object defining the orientation and ordering
#'   of the coordinate axes. If NULL, defaults to standard neurological convention
#'   (Left-Posterior-Inferior for 3D).
#' @param trans A transformation matrix mapping voxel indices to world coordinates. If NULL,
#'   constructed from spacing and origin.
#'
#' @return A new \code{\linkS4class{NeuroSpace}} object
#'
#' @importFrom methods new
#'
#' @section Coordinate Systems:
#' NeuroSpace manages two coordinate systems:
#' \itemize{
#'   \item Voxel coordinates: Zero-based indices into the image grid
#'   \item World coordinates: Real-world coordinates (typically in millimeters)
#' }
#'
#' The transformation between these systems is defined by:
#' \itemize{
#'   \item Voxel spacing (physical size of voxels)
#'   \item Origin (world coordinates of first voxel)
#'   \item Axis orientation (how image axes map to anatomical directions)
#' }
#'
#' @section Validation:
#' The constructor performs extensive validation:
#' \itemize{
#'   \item All dimensions must be positive integers
#'   \item All spacing values must be positive
#'   \item Origin and spacing must have matching lengths
#'   \item Transformation matrix must be invertible
#' }
#'
#' @examples
#' # Create a standard 3D space (64x64x40 voxels, 2mm isotropic)
#' space_3d <- NeuroSpace(
#'   dim = c(64L, 64L, 40L),
#'   spacing = c(2, 2, 2),
#'   origin = c(-90, -126, -72)
#' )
#'
#' # Check properties
#' dim(space_3d)           # Image dimensions
#' spacing(space_3d)       # Voxel sizes
#' origin(space_3d)        # World-space origin
#'
#' # Create a 2D slice space
#' space_2d <- NeuroSpace(
#'   dim = c(128L, 128L),
#'   spacing = c(1.5, 1.5),
#'   origin = c(-96, -96)
#' )
#'
#' # Convert between coordinate systems
#' world_coords <- c(0, 0, 0)
#' vox_idx <- coord_to_index(space_3d, world_coords)
#' back_to_world <- index_to_coord(space_3d, vox_idx)
#'
#' @seealso
#' \code{\linkS4class{AxisSet}} for axis orientation specification,
#' \code{\link{coord_to_index}} for coordinate conversion,
#' \code{\link{index_to_coord}} for inverse coordinate conversion,
#' \code{\linkS4class{NeuroObj}} for objects using NeuroSpace
#'
#' @references
#' For details on neuroimaging coordinate systems:
#' \itemize{
#'   \item Brett, M., Johnsrude, I. S., & Owen, A. M. (2002).
#'     The problem of functional localization in the human brain.
#'     Nature Reviews Neuroscience, 3(3), 243-249.
#'   \item Evans, A. C., et al. (1993). 3D statistical neuroanatomical models
#'     from 305 MRI volumes. Nuclear Science Symposium and Medical Imaging Conference.
#' }
#'
#' @export
NeuroSpace <- function(dim, spacing = NULL, origin = NULL, axes = NULL, trans = NULL) {
  if (!is.numeric(dim) || !all(dim == as.integer(dim)) || !all(dim > 0)) {
    cli::cli_abort("{.arg dim} must be a vector of positive integers.")
  }
  dim <- as.integer(dim)

  # Set defaults for spacing and origin
  if (is.null(spacing)) {
    spacing <- rep(1, min(length(dim), 3))
  }
  if (is.null(origin)) {
    origin <- rep(0, min(length(dim), 3))
  }

  if (!is.numeric(spacing) || !is.numeric(origin)) {
    cli::cli_abort("{.arg spacing} and {.arg origin} must be numeric vectors.")
  }
  if (length(origin) != length(spacing)) {
    cli::cli_abort("{.arg origin} and {.arg spacing} must have the same length.")
  }
  if (!all(spacing > 0)) {
    cli::cli_abort("All {.arg spacing} values must be positive.")
  }

  # Create transformation matrix if not provided
  if (is.null(trans)) {
    D <- min(length(dim), 3)
    trans <- diag(c(spacing, 1))
    trans[1:D, D+1] <- origin
  } else {
    if (!is.matrix(trans) || nrow(trans) != ncol(trans)) {
      cli::cli_abort("{.arg trans} must be a square matrix.")
    }
    # Derive spacing and origin from the affine when trans is provided
    D <- nrow(trans) - 1L
    spacing <- sqrt(colSums(trans[seq_len(D), seq_len(D), drop = FALSE]^2))
    origin <- trans[seq_len(D), D + 1L]
  }

  # Ensure matrix is invertible
  tryCatch({
    inverse <- solve(trans)
  }, error = function(e) {
    stop("transformation matrix must be invertible")
  })


  # Round to float32-level precision to avoid numerical issues
  # while preserving NIfTI-compatible accuracy (float32 ~ 7 digits)
  trans <- signif(trans, 7)
  inverse <- signif(inverse, 7)

  # Set up axes
  if (is.null(axes)) {
    if (length(dim) >= 3) {
      axes <- .nearestAnatomy(trans)
    } else if (length(dim) == 2) {
      axes <- AxisSet2D(LEFT_RIGHT, POST_ANT)
    } else {
      stop("unsupported number of dimensions")
    }
  }

  # Create object
  new("NeuroSpace",
      dim = dim,
      origin = signif(origin, 6),
      spacing = signif(spacing, 6),
      axes = axes,
      trans = trans,
      inverse = inverse)
}



#' add dimension to \code{\linkS4class{NeuroSpace}}
#' @param x The NeuroSpace object
#' @param n Numeric value specifying the size of the new dimension
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
            if (length(D) < 2) {
              cli::cli_abort("Cannot drop dimension from space with less than 2 dimensions.")
            }

            Dind <- seq_len(length(D))[-dimnum]
            if (ndim(x) > 3) {
              NeuroSpace(D[Dind], origin=origin(x)[Dind], spacing=spacing(x)[Dind],
                         axes=axes(x), trans=trans(x))
            } else {

              # For typical 2D/3D spaces, keep the rows/cols
              # corresponding to the retained spatial axes plus
              # the homogeneous coordinate, ensuring a square matrix.
              tx <- x@trans
              k  <- ndim(x)
              keep <- c(Dind, k + 1L)
              tx2 <- tx[keep, keep, drop=FALSE]

              NeuroSpace(D[Dind],
                         origin  = origin(x)[Dind],
                         spacing = spacing(x)[Dind],
                         axes    = drop_dim(axes(x), dimnum),
                         trans   = tx2)
            }

          })

#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x = "NeuroSpace", dimnum="missing"),
		def=function(x) {

			D <- dim(x)
			stopifnot(length(D) >= 2)
			Dind <- 1:(length(D)-1)


      if (ndim(x) > 3) {
			  NeuroSpace(D[Dind], origin=origin(x)[Dind], spacing=spacing(x)[Dind], axes=axes(x), trans=trans(x))
      } else {
        # Drop last spatial dimension from the affine:
        # keep rows/cols for retained axes + homogeneous coordinate
        tx <- trans(x)
        k  <- ndim(x)
        keep <- c(Dind, k + 1L)
        tx2 <- tx[keep, keep, drop = FALSE]
        NeuroSpace(D[Dind], axes=drop_dim(axes(x)), trans=tx2)
      }
		})


#' @export
#' @rdname dim-methods
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


#' Get dimension size along a specific axis
#' @param x The NeuroSpace object
#' @param axis The NamedAxis to query
#' @export
#' @rdname dim_of-methods
setMethod(f="dim_of", signature=signature(x = "NeuroSpace", axis="NamedAxis"),
            function(x, axis) {
  dir <- abs(axis@direction)

  dnum <- which_dim(x,axis)
  dim(x)[dnum]
})

#' Find which dimension corresponds to a given axis
#' @param x The NeuroSpace object
#' @param axis The NamedAxis to find
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
			} else if (ncol(coords) == 3) {
			  if (length(dx) < 3) {
			    cli::cli_abort("Space must have at least 3 dimensions for 3-column coordinate matrix, not {length(dx)}D.")
			  }
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
		    if (length(coords) != 2) {
		      cli::cli_abort("{.arg coords} must have length 2 for a 2D space, not {length(coords)}.")
		    }
		    if (coords[1] < 1 || coords[1] > dx[1] || coords[2] < 1 || coords[2] > dx[2]) {
		      cli::cli_abort("{.arg coords} ({.val {coords}}) are out of bounds for dimensions ({.val {dx}}).")
		    }
		    (coords[2]-1)*dx[1] + coords[1]
		  } else {
		    if (length(coords) != 3) {
		      cli::cli_abort("{.arg coords} must have length 3 for a 3D space, not {length(coords)}.")
		    }
			  .gridToIndex3D(dim(x), matrix(coords, nrow=1, byrow=TRUE))
		  }
		}
)



#' @export
#' @rdname coord_to_index-methods
setMethod(f="coord_to_index", signature=signature(x="NeuroVol", coords="matrix"),
          def=function(x, coords) {
            if (ncol(coords) != 3) {
              cli::cli_abort("Coordinate matrix must have 3 columns, not {ncol(coords)}.")
            }
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
            if (ncol(coords) != 3) {
              cli::cli_abort("Coordinate matrix must have 3 columns, not {ncol(coords)}.")
            }
            array.dim <- dim(x)
            .gridToIndex3D(dim(x), coords)
          })


#' @export
#' @rdname grid_to_index-methods
setMethod(f="grid_to_index", signature=signature(x="NeuroVol", coords="numeric"),
          def=function(x, coords) {
            if (length(coords) != 3) {
              cli::cli_abort("{.arg coords} must have length 3, not {length(coords)}.")
            }
            array.dim <- dim(x)
            .gridToIndex3D(dim(x), matrix(coords, nrow=1, byrow=TRUE))
          })


#' @export
#' @rdname reorient-methods
setMethod(f="reorient", signature=signature(x = "NeuroSpace", orient="character"),
          def=function(x, orient) {

            stopifnot(length(orient) == 3)
            anat <- findAnatomy3D(orient[1], orient[2], orient[3])
            pmat_new <- perm_mat(anat)


            tx <- t(pmat_new) %*% trans(x)[1:ndim(x),]
            tx <- rbind(tx,c(rep(0, ndim(x)),1))
            #itx <- zapsmall(MASS::ginv(tx))

            NeuroSpace(dim(x), spacing=spacing(x), axes=anat, trans=tx,
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
#' @rdname origin-methods
setMethod(f="origin", signature=signature(x = "NeuroVec"), def=function(x) space(x)@origin)


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


#' @export
#' @rdname show-methods
setMethod("show", "NeuroSpace", function(object) {
  d <- dim(object)
  show_header("NeuroSpace", paste0(ndim(object), "D"))
  show_rule("Geometry")
  show_field("Dimensions", paste(d, collapse = " x "))
  show_field("Spacing", paste(round(spacing(object), 3), collapse = " x "), " mm")
  show_field("Origin", paste(round(origin(object), 2), collapse = ", "))
  if (ndim(object) >= 3) {
    show_field("Orientation", safe_axcodes(object))
  }
  show_field("Voxels", format(prod(d), big.mark = ","))
})
