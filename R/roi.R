#' @include all_class.R
{}
#' @include all_generic.R
{}


#' Create an instance of class \code{\linkS4class{ROIVol}}
#'
#' This function constructs an instance of the ROIVol class, which represents
#' a region of interest (ROI) in a 3D volume. The class stores the
#' NeuroSpace object, voxel coordinates, and data values for the ROI.
#'
#' @param vspace An instance of class \code{NeuroSpace} with three dimensions,
#'   which represents the dimensions and voxel spacing of the 3D volume.
#' @param coords A 3-column matrix of voxel coordinates for the region of interest.
#' @param data The data values associated with the region of interest,
#'   provided as a numeric vector. By default, it is a vector of ones with a length equal
#'   to the number of rows in the `coords` matrix.
#' @return An instance of class \code{ROIVol}, containing the NeuroSpace object,
#'   voxel coordinates, and data values for the region of interest.
#' @examples
#' # Create a NeuroSpace object
#' vspace <- NeuroSpace(dim = c(5, 5, 5), spacing = c(1, 1, 1))
#'
#' # Define voxel coordinates for the ROI
#' coords <- matrix(c(1, 2, 3, 2, 2, 2, 3, 3, 3), ncol = 3)
#'
#' # Create a ROIVol object
#' roi_vol <- ROIVol(vspace, coords)
#' @export
ROIVol <- function(vspace, coords, data=rep(1, nrow(coords))) {
  new("ROIVol", space=vspace, coords=coords, as.vector(data))
}

#' Create an instance of class \code{\linkS4class{ROIVec}}
#'
#' This function constructs an instance of the ROIVec class, which represents
#' a region of interest (ROI) in a 4D volume. The class stores the
#' NeuroSpace object, voxel coordinates, and data values for the ROI.
#'
#' @param vspace An instance of class \code{NeuroSpace} with four dimensions,
#'   which represents the dimensions, voxel spacing, and time points of the 4D volume.
#' @param coords A 3-column matrix of voxel coordinates for the region of interest.
#' @param data The matrix of data values associated with the region of interest,
#'   with each row representing a voxel and each column representing a time point.
#'   By default, it is a matrix with a number of rows equal to the number of rows
#'   in the `coords` matrix and a single column filled with ones.
#' @return An instance of class \code{ROIVec}, containing the NeuroSpace object,
#'   voxel coordinates, and data values for the region of interest.
#' @examples
#' # Create a NeuroSpace object
#' vspace <- NeuroSpace(dim = c(5, 5, 5, 10), spacing = c(1, 1, 1))
#'
#' # Define voxel coordinates for the ROI
#' coords <- matrix(c(1, 2, 3, 2, 2, 2, 3, 3, 3), ncol = 3)
#'
#' # Create a data matrix for the ROI
#' data <- matrix(rnorm(30), nrow = 10, ncol = 3)
#'
#' # Create a ROIVec object
#' roi_vec <- ROIVec(vspace, coords, data)
#' @rdname ROIVec
#' @export
ROIVec <- function(vspace, coords, data=rep(nrow(coords),1)) {
  new("ROIVec", space=vspace, coords=coords, data)
}

#' convert a \code{\linkS4class{ROIVec}} to a matrix
#'
#' @rdname as.matrix-methods
#' @export
setMethod(f="as.matrix", signature=signature(x = "ROIVec"), def=function(x) {
  as(x, "matrix")
})


#' @keywords internal
#' @noRd
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

  if (all(map_int(coords, length) == 0)) {
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


#' @keywords internal
#' @noRd
.makeCubicGrid <- function(bvol, centroid, surround) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)

  coords <- lapply(centroid, function(x) { round(seq(x-surround, x+surround)) })
  coords <- lapply(1:3, function(i) {
    x <- coords[[i]]
    x[x > 0 & x <= vdim[i]]
  })

  if (all(map_int(coords, length) == 0)) {
    stop(paste("invalid cube for centroid", centroid, " with surround", surround, ": volume is zero"))
  }

  grid <- as.matrix(expand.grid(x=coords[[1]],y=coords[[2]],z=coords[[3]]))
}



#' Create a square region of interest
#'
#' This function creates a square region of interest (ROI) in a 3D volume, where the z-dimension is fixed
#' at one voxel coordinate. The ROI is defined within a given NeuroVol or NeuroSpace instance.
#'
#' @param bvol A \code{NeuroVol} or \code{NeuroSpace} instance representing the 3D volume or space.
#' @param centroid A numeric vector of length 3, representing the center of the square ROI in voxel coordinates.
#' @param surround A non-negative integer specifying the number of voxels on either side of the central voxel.
#' @param fill An optional value or values to assign to the data slot of the resulting ROI. If not provided, no data will be assigned.
#' @param nonzero A logical value indicating whether to keep only nonzero elements from \code{bvol}.
#'   If \code{bvol} is a \code{NeuroSpace} instance, this argument is ignored.
#' @param fixdim A logical value indicating whether the fixed dimension is the third, or z, dimension. Default is TRUE.
#' @return An instance of class \code{ROIVol} representing the square ROI.
#' @examples
#' sp1 <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
#' square <- square_roi(sp1, c(5, 5, 5), 1)
#' vox <- coords(square)
#' ## a 3 X 3 X 1 grid
#' nrow(vox) == 9
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
    vals != 0
  } else {
    TRUE
  }

  grid <- grid[keep,,drop=FALSE]

  center_index <- which(colSums(apply(grid, 1, "==", centroid)) == 3)
  parent_index <- grid_to_index(bvol, grid[center_index,])
  ### add central voxel
  new("ROIVolWindow", space=space(bvol), coords = grid, center_index=center_index, parent_index=parent_index,vals[keep])

}


#' Create A Cuboid Region of Interest
#'
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroid the center of the cube in \emph{voxel} coordinates
#' @param surround the number of voxels on either side of the central voxel. A \code{vector} of length 3.
#' @param fill optional value(s) to assign to data slot.
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{NeuroSpace} then this argument is ignored.
#' @return an instance of class \code{ROIVol}
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
#'  cube <- cuboid_roi(sp1, c(5,5,5), 3)
#'  vox <- coords(cube)
#'  cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=5)
#'
#'
#' @export
cuboid_roi <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }

  if (length(centroid) != 3) {
    stop("cuboid_roi: centroid must have length of 3 (x,y,z coordinates)")
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
    vals != 0
  } else {
    TRUE
  }

  grid <- grid[keep,,drop=FALSE]
  center_index <- which(colSums(apply(grid, 1, "==", centroid)) == 3)
  parent_index <- grid_to_index(bvol, grid[center_index,])

  new("ROIVolWindow", vals[keep], space=space(bvol), coords = grid, center_index=center_index, parent_index=parent_index)

}

#' @importFrom dbscan frNN
#' @keywords internal
#' @noRd
make_spherical_grid <- function(bvol, centroid, radius, use_cpp=TRUE) {

  vspacing <- spacing(bvol)

  if (radius < min(vspacing)) {
    stop("'radius' is too small; must be greater than at least one voxel dimension in image")
  }

  vdim <- dim(bvol)
  centroid <- as.integer(centroid)

  out <- if (use_cpp) {
    local_sphere(centroid[1], centroid[2], centroid[3], radius, vspacing, vdim)
  } else {
    deltas <- map_dbl(vspacing, function(x) round(radius/x))

    cube <- as.matrix(expand.grid(
      seq(centroid[1] - round(radius/vspacing[1]), centroid[1] + round(radius/vspacing[1])),
      seq(centroid[2] - round(radius/vspacing[2]), centroid[2] + round(radius/vspacing[2])),
      seq(centroid[3] - round(radius/vspacing[3]), centroid[3] + round(radius/vspacing[3]))))

    rs <- rowSums(sapply(1:ncol(cube), function(i) cube[,i] > 0 & cube[,i] <= vdim[i]))
    keep <- which(rs == 3)
    cube <- cube[keep,]

    coords <- t(t(cube) * vspacing)


    #res <- rflann::RadiusSearch(matrix(centroid * vspacing, ncol=3), coords, radius=radius^2,
    #                          max_neighbour=nrow(cube), build="kdtree", cores=0, checks=1)

    res <- dbscan::frNN(coords, eps=radius, query=matrix(centroid * vspacing, ncol=3))

    cube[res$id[[1]],,drop=FALSE]
  }

}

# masked_roi <- function(mask, vox, vox_offset) {
#   out <- t(vox + t(vox_offset))
#   d <- dim(mask)
#   keep <- (out[,1] > 0 & out[,1] < d[1]) & (out[,2] > 0 & out[,2] < d[2]) & (out[,3] > 0 & out[,3] < d[3])
#   out <- out[keep,]
#   vals <- mask[out]
#   out <- out[vals > 0,]
#   ROIVol(space(mask), out)
# }

#' Create a Spherical Region of Interest
#'
#' @description Creates a Spherical ROI based on a centroid.
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroid the center of the sphere in positive-coordinate (i,j,k) voxel space.
#' @param radius the radius in real units (e.g. millimeters) of the spherical ROI
#' @param fill optional value(s) to store as data
#' @param nonzero if \code{TRUE}, keep only nonzero elements from \code{bvol}
#' @param use_cpp whether to use compiled c++ code
#' @return an instance of class \code{ROIVol}
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))
#'  # create an ROI centered around the integer-valued positive voxel coordinate: i=5, j=5, k=5
#'  cube <- spherical_roi(sp1, c(5,5,5), 3.5)
#'  vox <- coords(cube)
#'  cds <- coords(cube, real=TRUE)
#'  ## fill in ROI with value of 6
#'  cube1 <- spherical_roi(sp1, c(5,5,5), 3.5, fill=6)
#'  all(cube1 == 6)
#'
#'  # create an ROI centered around the real-valued coordinates: x=5, y=5, z=5
#'  vox <- coord_to_grid(sp1, c(5, 5, 5))
#'  cube <- spherical_roi(sp1, vox, 3.5)
#' @export
spherical_roi <- function (bvol, centroid, radius, fill=NULL, nonzero=FALSE, use_cpp=TRUE) {
  if (is.matrix(centroid)) {
    assertthat::assert_that(ncol(centroid == 3) & nrow(centroid) == 1)
    centroid <- drop(centroid)
  }

  vdim <- dim(bvol)

  assertthat::assert_that(length(centroid) == 3)
  assertthat::assert_that(centroid[1] <= vdim[1] && centroid[2] <= vdim[2] && centroid[3] <= vdim[3])
  assertthat::assert_that(all(centroid > 0))

  if (is.null(fill) && is(bvol, "NeuroSpace")) {
    fill = 1
  }

  bspace <- space(bvol)
  vspacing <- spacing(bvol)

  centroid <- as.integer(centroid)
  grid <- make_spherical_grid(bvol, centroid, radius, use_cpp=use_cpp)

  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    as.numeric(bvol[grid])
  }

  if (nonzero) {
    keep <- vals != 0
    grid <- grid[keep, ,drop=FALSE]
    center_index <- which(colSums(apply(grid, 1, "==", centroid)) == 3)
    parent_index <- grid_to_index(bvol, grid[center_index,])
    new("ROIVolWindow", vals[keep], space=bspace, coords = grid, center_index=center_index, parent_index=parent_index)
  } else {
    center_index <- which(colSums(apply(grid, 1, "==", centroid)) == 3)
    parent_index <- grid_to_index(bvol, grid[center_index,])
    new("ROIVolWindow", vals, space=bspace, coords = grid,center_index=center_index, parent_index=parent_index)
  }

}

# spherical_basis <- function(bvol, coord, kernel, weight=1) {
#   ## convert coordinate from MNI space to voxel space
#   grid.loc <- coord_to_grid(bvol, coord)
#
#   ## shift kernel so that it is centered around 'grid.loc'
#   voxmat <- floor(voxels(kernel, centerVoxel=grid.loc))
#   indices <- gridToIndex(template, voxmat)
#   neuroim:::SparseNeuroVol(kernel@weights * weight, template, indices=indices)
# }

#' @keywords internal
#' @noRd
.resample <- function(x, ...) x[sample.int(length(x), ...)]

#' @keywords internal
#' @noRd
roi_vector_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_vector_matrix", "matrix"))

}

#' @keywords internal
#' @noRd
roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))

}



#' Coerce ROIVec to matrix
#'
#' This function provides a method to coerce an object of class \code{ROIVec} to a \code{matrix}.
#'
#' @name as
#' @param from An object of class \code{ROIVec} to be coerced to a \code{matrix}.
#' @return A \code{matrix} obtained by coercing the \code{ROIVec} object.
#' @export
setAs(from="ROIVec", to="matrix", function(from) {
  ind <- indices(from)
  roi_vector_matrix(from@.Data, refspace=from@space, indices=ind,
                    coords=index_to_coord(drop_dim(from@space),
                                          as.numeric(ind)))

})

#' Coerce ROIVol to DenseNeuroVol
#'
#' This function provides a method to coerce an object of class \code{ROIVol} to a \code{DenseNeuroVol}.
#'
#' @name as
#' @param from An object of class \code{ROIVol} to be coerced to a \code{DenseNeuroVol}.
#' @return A \code{DenseNeuroVol} object obtained by coercing the \code{ROIVol} object.
setAs(from="ROIVol", to="DenseNeuroVol", function(from) {
  NeuroVol(values(from), space(from), indices=indices(from))
  #dat <- array(0, dim(from@space))
  #dat[coords(from)] <- from@data
  #ovol <- DenseNeuroVol(dat, from@space, from@source)
})


#' Coerce ROIVol to DenseNeuroVol using as.dense method
#'
#' This function provides a method to coerce an object of class \code{ROIVol} to a \code{DenseNeuroVol} using the \code{as.dense} method.
#'
#' @rdname as.dense-methods
#' @param x An object of class \code{ROIVol} to be coerced to a \code{DenseNeuroVol}.
#' @return A \code{DenseNeuroVol} object obtained by coercing the \code{ROIVol} object.
setMethod("as.dense", signature(x="ROIVol"),
          function(x) {
            as(x, "DenseNeuroVol")
            #NeuroVol(values(x), space(x), indices=indices(x))
})

#' @export
#' @rdname centroid-methods
setMethod(f="centroid", signature=signature(x = "ROICoords"),
          def=function(x)  {
            cds = coords(x, real=TRUE)
            colMeans(cds)
          })



#' @rdname values-methods
#' @export
setMethod("values", signature(x="ROIVol"),
          function(x, ...) {
             x@.Data
          })

#' @rdname values-methods
#' @export
setMethod("values", signature(x="ROIVec"),
          function(x, ...) {
            x@.Data
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="ROIVec", subset="missing"),
          function(x) {
            ind <- 1:nrow(x@coords)
            f <- function(i) x@.Data[,i]
            #lis <- map(ind, function(i) f)
            deflist::deflist(f, length(ind))
          })


#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="matrix", subset="missing"),
          function(x) {
            ind <- 1:ncol(x)
            f <- function(i) x[,i]
            #lis <- map(ind, function(i) f)
            deflist::deflist(f, length(ind))
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="ROIVec", subset="integer"),
          function(x, subset) {
            ind <- (1:nrow(x@coords))[subset]
            f <- function(i) x@.Data[,ind[i]]
            #lis <- map(ind, function(i) f)
            deflist::deflist(f, length(ind))
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="matrix", subset="integer"),
          function(x, subset) {
            ind <- (1:ncol(x))[subset]
            f <- function(i) x[,ind[i]]
            #lis <- map(ind, function(i) f)
            deflist::deflist(f, length(ind))
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="matrix", subset="numeric"),
          function(x, subset) {
            callGeneric(x,subset)
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="ROIVec", subset="numeric"),
          function(x, subset) {
            callGeneric(x, as.integer(subset))
          })

#' @rdname vectors-methods
#' @export
setMethod("vectors", signature(x="ROIVec", subset="logical"),
          function(x, subset) {
            callGeneric(x, as.integer(which(subset)))
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
            .gridToIndex(dim(x@space)[1:3], x@coords)
            #grid_to_index(x@space, x@coords)
          })


#' @export
#' @param real if \code{TRUE}, return coordinates in real world units
#' @rdname coords-methods
setMethod(f="coords", signature=signature(x="ROICoords"),
          function(x, real=FALSE) {
            if (real) {
              input <- t(cbind(x@coords-.5, rep(1, nrow(x@coords))))
              ret <- t(trans(x@space) %*% input)
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
            ROIVol(x@space, x@coords[i,,drop=FALSE], x@.Data[i])
          })

#' @rdname vol_subset-methods
#' @aliases [,ROIVol,logical,missing,ANY-method
setMethod("[", signature=signature(x="ROIVol", i="logical", j="missing", drop="ANY"),
          function(x,i,j,drop) {
            ROIVol(x@space, x@coords[i,,drop=FALSE], x@.Data[i])
          })


#' show an \code{\linkS4class{ROIVol}}
#' @param object the object
#' @export
setMethod("show", signature=signature(object = "ROIVol"),
		  function (object) {
			  cat("\n\nROIVol", "\n")
			  cat("  Size:           ", length(object), "\n")
			  cat("  Parent Dim:     ", dim(object@space), "\n")
			  cat("  Num Data Cols:  ", 1, "\n" )
			  cat("  Voxel Cen. Mass:", colMeans(coords(object)), "\n")
		  })


#' show an \code{\linkS4class{ROIVec}}
#' @param object the object
#' @export
setMethod("show", signature=signature(object = "ROIVec"),
          function (object) {
            cat("\n\nROIVec", "\n")
            cat("  ncol:           ", ncol(object), "\n")
            cat("  nrow:           ", nrow(object), "\n")
            cat("  Parent Dim:     ", dim(object@space), "\n")
            cat("  Voxel Cen. Mass:", colMeans(coords(object)), "\n")
          })



#' Create a Kernel object from a function of distance from kernel center
#'
#' This function creates a Kernel object using a kernel function (\code{FUN}) that takes the distance from the center of the kernel as its first argument.
#'
#' @param kerndim A numeric vector representing the dimensions in voxels of the kernel.
#' @param vdim A numeric vector representing the dimensions of the voxels in real units.
#' @param FUN The kernel function taking its first argument representing the distance from the center of the kernel (default: \code{dnorm}).
#' @param ... Additional parameters to the kernel function, \code{FUN}.
#' @importFrom stats dnorm
#' @return A Kernel object with the specified dimensions, voxel dimensions, and kernel function.
#' @examples
#' kdim <- c(3, 3, 3)
#' vdim <- c(1, 1, 1)
#' k <- Kernel(kerndim = kdim, vdim = vdim, FUN = dnorm, sd = 1)
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
  hwidth <- map_dbl(kerndim, function(d) ceiling(d/2 -1))

  ## note, if a kernel dim is even, this will force it to be odd numbered
  grid.vec <- map(hwidth, function(sv) seq(-sv, sv))

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

#' @export
#' @rdname embed_kernel-methods
#' @param weight multiply kernel weights by this value
setMethod("embed_kernel", signature=signature(x="Kernel", sp="NeuroSpace", center_voxel="numeric"),
          function(x,  sp, center_voxel, weight=1) {
            vox <- floor(voxels(x, center_voxel))
            indices <- grid_to_index(sp, vox)
            SparseNeuroVol(x@weights * weight, sp, indices=indices)
          })



#' @param center_voxel the absolute location of the center of the voxel, default is (0,0,0)
#' @rdname voxels-methods
#' @export
setMethod(f="voxels", signature=signature(x="Kernel"),
          function(x, center_voxel=NULL) {
            if (is.null(center_voxel)) {
              x@voxels
            } else {
              sweep(x@voxels, 2, center_voxel, "+")
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


