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
#' @return A new \code{NeuroSpace} object with the following slots:
#' \itemize{
#'   \item \code{dim}: Integer vector of image dimensions
#'   \item \code{spacing}: Numeric vector of voxel sizes
#'   \item \code{origin}: Numeric vector of world-space origin
#'   \item \code{axes}: AxisSet object defining orientation
#'   \item \code{trans}: Matrix for voxel-to-world transformation
#'   \item \code{inverse}: Matrix for world-to-voxel transformation
#' }
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
#' \itemize{
#'   \item \code{\link{AxisSet}} for axis orientation specification
#'   \item \code{\link{coord_to_index}} for coordinate conversion
#'   \item \code{\link{index_to_coord}} for inverse coordinate conversion
#'   \item \code{\link{NeuroObj}} for objects using NeuroSpace
#' }
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
#' @importFrom assertthat assert_that
NeuroSpace <- function(dim, spacing = NULL, origin = NULL, axes = NULL, trans = NULL) {
  # Input validation
  if (!is.numeric(dim) || any(dim != as.integer(dim)) || any(dim <= 0)) {
    stop("'dim' must be a vector of positive integers")
  }
  dim <- as.integer(dim)
  
  # Set defaults for spacing and origin
  if (is.null(spacing)) {
    spacing <- rep(1, min(length(dim), 3))
  }
  if (is.null(origin)) {
    origin <- rep(0, min(length(dim), 3))
  }
  
  # Validate spacing and origin
  if (!is.numeric(spacing) || !is.numeric(origin)) {
    stop("'spacing' and 'origin' must be numeric vectors")
  }
  if (length(origin) != length(spacing)) {
    stop("'origin' and 'spacing' must have the same length")
  }
  if (any(spacing <= 0)) {
    stop("all 'spacing' values must be positive")
  }
  
  # Create transformation matrix if not provided
  if (is.null(trans)) {
    D <- min(length(dim), 3)
    trans <- diag(c(spacing, 1))
    trans[1:D, D+1] <- origin
  } else if (!is.matrix(trans) || nrow(trans) != ncol(trans)) {
    stop("'trans' must be a square matrix")
  }
  
  # Ensure matrix is invertible
  tryCatch({
    inverse <- solve(trans)
  }, error = function(e) {
    stop("transformation matrix must be invertible")
  })
  
  # Round to avoid numerical issues
  trans <- signif(trans, 6)
  inverse <- signif(inverse, 6)
  
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

#' Display NeuroSpace Object Details
#'
#' @description
#' Provides a beautifully formatted display of a NeuroSpace object's properties,
#' including dimensions, spacing, origin, and transformation information.
#' Uses color coding and unicode characters for enhanced readability.
#'
#' @param object A NeuroSpace object to display
#'
#' @importFrom crayon bold blue green red yellow silver white bgBlue
#' @importFrom utils object.size
#' @export
setMethod(f="show", 
          signature=signature("NeuroSpace"),
          def=function(object) {
            # Helper function for memory size formatting
            format_bytes <- function(bytes) {
              units <- c('B', 'KB', 'MB', 'GB', 'TB')
              i <- floor(log2(bytes) / 10)
              sprintf("%.1f %s", bytes / 2^(10 * i), units[i + 1])
            }
            
            # Helper for matrix formatting
            format_matrix <- function(mat, digits=3) {
              formatted <- apply(mat, 1, function(row) {
                paste(sprintf(paste0("%", digits + 4, ".", digits, "f"), row), collapse="  ")
              })
              paste(formatted, collapse="\n")
            }
            
            # Calculate memory footprint
            mem_size <- format_bytes(utils::object.size(object))
            
            # Header
            cat("\n")
            cat(bgBlue(white(bold(" NeuroSpace Object "))))
            cat("\n")
            
            # Dimension Information
            cat("\n", bold(yellow("▶ Dimensions")), "\n")
            dim_str <- paste(object@dim, collapse=" × ")
            cat("  ", silver("Grid Size:"), " ", green(dim_str), "\n", sep="")
            cat("  ", silver("Memory:"), "   ", green(mem_size), "\n", sep="")
            
            # Spatial Properties
            cat("\n", bold(yellow("▶ Spatial Properties")), "\n")
            spacing_str <- paste(sprintf("%.2f", object@spacing), collapse=" × ")
            origin_str <- paste(sprintf("%.2f", object@origin), collapse=" × ")
            cat("  ", silver("Spacing:"), "   ", blue(spacing_str), " ", silver("mm"), "\n", sep="")
            cat("  ", silver("Origin:"), "    ", blue(origin_str), " ", silver("mm"), "\n", sep="")
            
            # Anatomical Orientation
            cat("\n", bold(yellow("▶ Anatomical Orientation")), "\n")
            if (length(object@dim) >= 3) {
              orientations <- c(
                paste0("X: ", green(object@axes@i@axis)),
                paste0("Y: ", green(object@axes@j@axis)),
                paste0("Z: ", green(object@axes@k@axis))
              )
            } else {
              orientations <- c(
                paste0("X: ", green(object@axes@i@axis)),
                paste0("Y: ", green(object@axes@j@axis))
              )
            }
            cat(paste0("  ", paste(orientations, collapse="  |  ")), "\n")
            
            # World Transformation
            cat("\n", bold(yellow("▶ World Transformation")), "\n")
            cat(silver("  Forward (Voxel → World):"), "\n")
            cat(blue(paste0("    ", format_matrix(object@trans))), "\n")
            cat(silver("  Inverse (World → Voxel):"), "\n")
            cat(blue(paste0("    ", format_matrix(object@inverse))), "\n")
            
            # Bounding Box (in world coordinates)
            cat("\n", bold(yellow("▶ Bounding Box")), "\n")
            ndim <- length(object@dim)
            if (ndim == 2) {
              corners <- matrix(c(0, 0, 
                                object@dim[1]-1, 0,
                                0, object@dim[2]-1,
                                object@dim[1]-1, object@dim[2]-1), 
                              nrow=4, byrow=TRUE)
              world_corners <- t(object@trans[1:2, 1:2, drop=FALSE] %*% t(corners) + 
                                 matrix(object@trans[1:2, 3], nrow=2, ncol=4))
              min_corner <- apply(world_corners, 2, min)
              max_corner <- apply(world_corners, 2, max)
              cat("  ", silver("Min Corner:"), " ", 
                  green(paste(sprintf("%.1f", min_corner), collapse=", ")), 
                  " mm\n", sep="")
              cat("  ", silver("Max Corner:"), " ", 
                  green(paste(sprintf("%.1f", max_corner), collapse=", ")), 
                  " mm\n", sep="")
            } else {
              corners <- matrix(c(0, 0, 0, 
                                object@dim[1]-1, 0, 0,
                                0, object@dim[2]-1, 0,
                                object@dim[1]-1, object@dim[2]-1, 0,
                                0, 0, object@dim[3]-1,
                                object@dim[1]-1, 0, object@dim[3]-1,
                                0, object@dim[2]-1, object@dim[3]-1,
                                object@dim[1]-1, object@dim[2]-1, object@dim[3]-1), 
                              nrow=8, byrow=TRUE)
              world_corners <- t(object@trans[1:3, 1:3, drop=FALSE] %*% t(corners) + 
                                 matrix(object@trans[1:3, 4], nrow=3, ncol=8))
              min_corner <- apply(world_corners, 2, min)
              max_corner <- apply(world_corners, 2, max)
              cat("  ", silver("Min Corner:"), " ", 
                  green(paste(sprintf("%.1f", min_corner), collapse=", ")), 
                  " mm\n", sep="")
              cat("  ", silver("Max Corner:"), " ", 
                  green(paste(sprintf("%.1f", max_corner), collapse=", ")), 
                  " mm\n", sep="")
            }
            
            # Footer
            cat("\n", bgBlue(white(bold(paste(rep("═", 50), collapse="")))), "\n", sep="")
          })

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
