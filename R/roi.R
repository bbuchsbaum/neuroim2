#' @import iterators

#' @include all_class.R
{}
#' @include all_generic.R
{}


#' Create an instance of class \code{\linkS4class{ROIVol}}
#'
#' @param space an instance of class \code{NeuroSpace}
#' @param coords matrix of voxel coordinates
#' @param data the data values, numeric vector
#' @return an instance of class \code{ROIVol}
#' @rdname ROIVol
#' @export
ROIVol <- function(vspace, coords, data=rep(nrow(coords),1)) {
  new("ROIVol", space=vspace, coords=coords, data=as.vector(data))
}

#' Create an instance of class \code{\linkS4class{ROIVec}}
#'
#' @param vspace an instance of class \code{NeuroSpace}
#' @param coords matrix of voxel coordinates
#' @param data the \code{matrix} of data values
#' @return an instance of class \code{ROIVec}
#' @rdname ROIVec
#' @export
ROIVec <- function(vspace, coords, data=rep(nrow(coords),1)) {
  new("ROIVec", space=vspace, coords=coords, data=data)
}

#' convert a \code{ROIVec} to a matrix
#'
#' @rdname as.matrix-methods
#' @param x the object
#' @export
setMethod(f="as.matrix", signature=signature(x = "ROIVec"), def=function(x) {
  as(x, "matrix")
})



.makeSquareGrid <- function(bvol, centroid, surround, fixdim=3) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)

  dimnums <- seq(1,3)[-fixdim]

  coords <- lapply(centroid, function(x) { round(seq(x-surround, x+surround)) })
  coords <- lapply(dimnums, function(i) {
    x <- coords[[i]]
    x[x > 0 & x <= vdim[i]]
  })

  if (all(sapply(coords, length) == 0)) {
    stop(paste("invalid cube for centroid", centroid, " with surround", surround, ": volume is zero"))
  }

  if (fixdim == 3) {
    grid <- as.matrix(expand.grid(x=coords[[1]],y=coords[[2]],z=centroid[3]))
  } else if (fixdim == 2) {
    grid <- as.matrix(expand.grid(x=coords[[1]],y=centroid[2],z=coords[[2]]))
  } else if (fixdim == 1) {
    grid <- as.matrix(expand.grid(x=centroid[1],y=coords[[1]],z=coords[[2]]))
  }

  grid

}

.makeCubicGrid <- function(bvol, centroid, surround) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)

  coords <- lapply(centroid, function(x) { round(seq(x-surround, x+surround)) })
  coords <- lapply(1:3, function(i) {
    x <- coords[[i]]
    x[x > 0 & x <= vdim[i]]
  })

  if (all(sapply(coords, length) == 0)) {
    stop(paste("invalid cube for centroid", centroid, " with surround", surround, ": volume is zero"))
  }

  grid <- as.matrix(expand.grid(x=coords[[1]],y=coords[[2]],z=coords[[3]]))
}



#' Create a square region of interest where the z-dimension is fixed at one voxel coordinate.
#'
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance.
#' @param centroid the center of the cube in \emph{voxel} coordinates.
#' @param surround the number of voxels on either side of the central voxel.
#' @param fill optional value(s) to assign to data slot.
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{NeuroSpace} then this argument is ignored.
#' @param fixdim the fixed dimension is the third, or z, dimension.
#' @return an instance of class \code{ROIVol}.
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
#'  square <- square_roi(sp1, c(5,5,5), 1)
#'  vox <- coords(square)
#'  ## a 3 X 3 X 1 grid
#'  nrow(vox) == 9
#' @export
square_roi <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE, fixdim=3) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }

  if (length(centroid) != 3) {
    stop("square_roi: centroid must have length of 3 (x,y,z coordinates)")
  }

  if (surround < 0) {
    stop("'surround' argument cannot be negative")
  }

  if (is(bvol, "NeuroSpace") && is.null(fill)) {
    fill = 1
  }

  grid <- .makeSquareGrid(bvol,centroid,surround,fixdim=fixdim)

  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    as.numeric(bvol[grid])
  }

  keep <- if (nonzero) {
    which(vals != 0)
  } else {
    seq_along(vals)
  }

  ### add central voxel
  ROIVol(space(bvol), data = vals[keep], coords = grid[keep, ])

}


#' Create A Cuboid Region of Interest
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroid the center of the cube in \emph{voxel} coordinates
#' @param surround the number of voxels on either side of the central voxel. A \code{vector} of length 3.
#' @param fill optional value(s) to assign to data slot.
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{NeuroSpace} then this argument is ignored.
#' @return an instance of class \code{ROIVol}
#' @rdname cube_roi
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
#'  cube <- cube_roi(sp1, c(5,5,5), 3)
#'  vox <- coords(cube)
#'  cube2 <- cube_roi(sp1, c(5,5,5), 3, fill=5)
#'
#'
#' @export
cube_roi <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }

  if (length(centroid) != 3) {
    stop("cube_roi: centroid must have length of 3 (x,y,z coordinates)")
  }

  if (surround < 0) {
    stop("'surround' argument cannot be negative")
  }

  if (is(bvol, "NeuroSpace") && is.null(fill)) {
    fill = 1
  }

  grid <- .makeCubicGrid(bvol,centroid,surround)

  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    as.numeric(bvol[grid])
  }

  keep <- if (nonzero) {
    which(vals != 0)
  } else {
    seq_along(vals)
  }

  ### add central voxel
  ROIVol(space(bvol), data = vals[keep], coords = grid[keep, ])

}

#' @importFrom rflann RadiusSearch
.makeSphericalGrid <- function(bvol, centroid, radius) {

  vspacing <- spacing(bvol)

  if (radius < min(vspacing)) {
    stop("'radius' is too small; must be greater than at least one voxel dimension in image")
  }

  vdim <- dim(bvol)
  centroid <- as.integer(centroid)


  cube <- as.matrix(expand.grid(
    seq(centroid[1] - round(radius/vspacing[1]), centroid[1] + round(radius/vspacing[1])),
    seq(centroid[2] - round(radius/vspacing[2]), centroid[2] + round(radius/vspacing[2])),
    seq(centroid[3] - round(radius/vspacing[3]), centroid[3] + round(radius/vspacing[3]))))


  coords <- t(t(cube) * vspacing)

  res <- rflann::RadiusSearch(matrix(centroid * vspacing, ncol=3), coords, radius=radius^2, max_neighbour=nrow(cube))

  cube[res$indices[[1]],]

}



#' @title Create a Spherical Region of Interest
#'
#' @description Creates a Spherical ROI based on a Centroid.
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroid the center of the sphere in voxel space
#' @param radius the radius in real units (e.g. millimeters) of the spherical ROI
#' @param fill optional value(s) to store as data
#' @param nonzero if \code{TRUE}, keep only nonzero elements from \code{bvol}
#' @return an instance of class \code{ROIVol}
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))
#'  cube <- spherical_roi(sp1, c(5,5,5), 3.5)
#'  vox <- coords(cube)
#'  cds <- coords(cube, real=TRUE)
#'  ## fill in ROI with value of 6
#'  cube1 <- spherical_roi(sp1, c(5,5,5), 3.5, fill=6)
#'  all(cube1@data == 6)
#' @export
spherical_roi <- function (bvol, centroid, radius, fill=NULL, nonzero=FALSE) {
  if (is.matrix(centroid)) {
    assertthat::assert_that(ncol(centroid == 3) & nrow(centroid) == 1)
    centroid <- drop(centroid)
  }

  assertthat::assert_that(length(centroid) == 3)

  if (is.null(fill) && is(bvol, "NeuroSpace")) {
    fill = 1
  }

  bspace <- space(bvol)
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)
  grid <- .makeSphericalGrid(bvol, centroid, radius)

  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    as.numeric(bvol[grid])
  }

  if (nonzero) {
    keep <- vals != 0
    ROIVol(bspace, data = vals[keep], coords = grid[keep, ,drop=FALSE])
  } else {
    ROIVol(bspace, data = vals, coords = grid)
  }

}

.resample <- function(x, ...) x[sample.int(length(x), ...)]

roi_vector_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_vector_matrix", "matrix"))

}

roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))

}



#' @name as
#' @rdname as-methods
setAs(from="ROIVec", to="matrix", function(from) {
  ind <- indices(from)
  roi_vector_matrix(from@data, refspace=from@space, indices=ind,
                    coords=index_to_coord(drop_dim(from@space),
                                          as.numeric(ind)))

})


#' @name as
#' @rdname as-methods
setAs(from="ROIVol", to="DenseNeuroVol", function(from) {
  dat <- array(0, dim(from@space))
  dat[coords(from)] <- from@data
  ovol <- DenseNeuroVol(dat, from@space, from@source)
})


#' @rdname values-methods
#' @export
setMethod("values", signature(x="ROIVol"),
          function(x, ...) {
             x@data
          })

#' @rdname values-methods
#' @export
setMethod("values", signature(x="ROIVec"),
          function(x, ...) {
            x@data
          })


#' @rdname indices-methods
#' @export
setMethod("indices", signature(x="ROIVol"),
          function(x) {
			      grid_to_index(x@space, x@coords)
		  })

#' @rdname indices-methods
#' @export
setMethod("indices", signature(x="ROIVec"),
          function(x) {
            grid_to_index(x@space, x@coords)
          })


#' @export
#' @param real if \code{TRUE}, return coordinates in real world units
#' @rdname coords-methods
setMethod(f="coords", signature=signature(x="ROIVol"),
          function(x, real=FALSE) {
            if (real) {
              input <- t(cbind(x@coords-.5, rep(1, nrow(x@coords))))
              ret <- t(trans(x) %*% input)
              ret[,1:3,drop=FALSE]
            } else {
              x@coords
            }
          })


#' @export
#' @rdname length-methods
setMethod(f="length", signature=signature(x="ROIVol"),
          function(x) {
            nrow(x@coords)
          })





#' subset an \code{ROIVol}
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param drop drop dimension
#' @rdname vol_subset-methods
#' @aliases [,ROIVol,numeric,missing,ANY-method
setMethod("[", signature=signature(x = "ROIVol", i = "numeric", j = "missing", drop = "ANY"),
          function (x, i, j, drop) {
            ROIVol(x@space, x@coords[i,,drop=FALSE], x@data[i])
          })

#' @rdname vol_subset-methods
#' @aliases [,ROIVol,logical,missing,ANY-method
setMethod("[", signature=signature(x="ROIVol", i="logical", j="missing", drop="ANY"),
          function(x,i,j,drop) {
            ROIVol(x@space, x@coords[i,,drop=FALSE], x@data[i])
          })


#' show an \code{\linkS4class{ROIVol}}
#' @param object the object
#' @export
setMethod("show", signature=signature(object = "ROIVol"),
		  function (object) {
			  cat("\n\n\tROIVol", "\n")
			  cat("\t size: ", length(object), "\n")
			  cat("\t parent dim:", dim(object), "\n")
			  cat("\t num data cols:", if (is.matrix(object@data)) ncol(object@data) else 1, "\n" )
			  cat("\t voxel center of mass: ", colMeans(coords(object)), "\n")
		  })








#' Create a Kernel object from a function of distance from kernel center
#'
#' @param kerndim the dimensions in voxels of the kernel
#' @param vdim the dimensions of the voxels in real units
#' @param FUN the kernel function taking as its first argument representing the distance from the center of the kernel
#' @param ... additional parameters to the kernel FUN
#' @importFrom stats dnorm
#' @export
Kernel <- function(kerndim, vdim, FUN=dnorm, ...) {
  if (length(kerndim) < 2) {
    stop("kernel dim length must be greater than 1")
  }

  .distance <- function(p1, p2) {
    diffs = (p1 - p2)
    sqrt(sum(diffs*diffs))
  }

  #kern <- array(0, kerndim)

  ## the half-width for each dimensions
  hwidth <- sapply(kerndim, function(d) ceiling(d/2 -1))

  ## note, if a kernel dim is even, this will force it to be odd numbered
  grid.vec <- lapply(hwidth, function(sv) seq(-sv, sv))

  # compute relative voxel locations (i.e. centered at 0,0,0)
  voxel.ind <- as.matrix(do.call("expand.grid", grid.vec))

  # fractional voxel locations so that the location of a voxel coordinate is centered within the voxel
  cvoxel.ind <- t(apply(voxel.ind, 1, function(vals) sign(vals)* ifelse(vals == 0, 0, abs(vals)-.5)))

  ## the coordinates ofthe voxels (i.e. after multiplying by pixel dims)
  coords <- t(apply(cvoxel.ind, 1, function(v) (v * vdim)))

  ## distance of coordinate from kernel center
  coord.dist <- apply(coords, 1, .distance, c(0,0,0))

  wts <- FUN(coord.dist, ...)
  wts <- wts/sum(wts)


  kern.weights <- wts

  new("Kernel", width=kerndim, weights=kern.weights, voxels=voxel.ind, coords=coords)

}



#' @param centerVoxel the absolute location of the center of the voxel, default is (0,0,0)
#' @rdname voxels-methods
#' @export
setMethod(f="voxels", signature=signature(x="Kernel"),
          function(x, centerVoxel=NULL) {
            if (is.null(centerVoxel)) {
              x@voxels
            } else {
              sweep(x@voxels, 2, centerVoxel, "+")
            }
          })



# GradientKernel <- function(direction=c("x", "y", "z")) {
#   direction <- match.arg(direction)
#   grid.vec <- lapply(1:3, function(sv) seq(-1, 1))
#
#   # compute relative voxel locations (i.e. centered at 0,0,0)
#   voxel.ind <- as.matrix(do.call("expand.grid", grid.vec))
#
#   # fractional voxel locations so that the location of a voxel coordinate is centered within the voxel
#   cvoxel.ind <- t(apply(voxel.ind, 1, function(vals) sign(vals)* ifelse(vals == 0, 0, abs(vals)-.5)))
#
#   ## the coordinates of the voxels (i.e. after multiplying by pixel dims)
#   coords <- t(apply(cvoxel.ind, 1, function(v) (v * vdim)))
#
#   if (direction == "x") {
#     gdim <- 1
#     odim <- 2:3
#   } else if (direction == "y") {
#     gdim <- 2
#     odim <- c(1,3)
#   } else {
#     gdim <- 3
#     odim <- c(1,2)
#   }
#
#
#   wts <- apply(coords, 1, function(r) {
#     if (r[gdim] == 0) {
#       0
#     } else if (r[gdim] < 0 && all(r[odim] == 0)) {
#       -2
#     } else if (r[gdim] > 0 && all(r[odim] == 0)) {
#       2
#     } else if (r[gdim] < 0) {
#       -1
#     } else {
#       1
#     }
#   })
#
#   new("Kernel", width=c(3,3,3), weights=wts, voxels=voxel.ind, coords=coords)
#
# }
#


