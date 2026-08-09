#' @include all_class.R
NULL
#' @include all_generic.R
NULL
#' @importFrom methods new
NULL


# ---------------------------------------------------------------------------
# Fast ROI construction
#
# `new("ROIVolWindow", ...)` costs ~450us and the cost is independent of ROI
# size: the class reaches a basic type ("numeric") through two levels of
# contains=, so the default initialize() walks that chain with callNextMethod()
# and runs validObject() at each level -- and validObject() recurses into the
# NeuroSpace slot and its nested AxisSet/NamedAxis tree. In an exhaustive
# searchlight that fixed overhead is ~70% of total runtime.
#
# Cloning a cached prototype and assigning the slots produces an object that is
# identical() to the new() result for ~1/11th of the cost. The class invariants
# that validObject() would have checked are asserted directly here instead, so
# nothing is silently skipped -- they are just checked once, cheaply, in the one
# place these objects are built.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.roi_proto_cache <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.roi_prototype <- function(class_name) {
  hit <- .roi_proto_cache[[class_name]]
  if (!is.null(hit)) {
    return(hit)
  }

  sp <- NeuroSpace(c(1L, 1L, 1L))
  proto <- switch(
    class_name,
    ROIVolWindow = new("ROIVolWindow", numeric(0), space = sp,
                       coords = matrix(0, nrow = 0, ncol = 3),
                       center_index = integer(0), parent_index = NA_integer_),
    cli::cli_abort("No ROI prototype for {.cls {class_name}}.")
  )

  assign(class_name, proto, envir = .roi_proto_cache)
  proto
}

#' Prepare a value for the .Data part the way new() would
#'
#' \code{new()} coerces whatever it is handed into the basic type the class
#' extends and drops the attributes -- including \code{names} -- that the data
#' part must not carry. Slot assignment does neither, so an ROI built by the
#' fast constructor from, say, \code{fill = TRUE} or \code{fill = c(a = 1)}
#' would otherwise differ from one built by \code{new()}. Integer storage is
#' left alone because \code{new()} preserves it too.
#'
#' @keywords internal
#' @noRd
.as_roi_data <- function(data) {
  if (is.null(data)) {
    stop("cannot use object of class \"NULL\" in new()")
  }
  tp <- typeof(data)
  # is.object() keeps classed values (notably factors, which are integer
  # underneath) out of the pass-through branch: as.vector() on a factor yields
  # its labels, whereas new() coerces through the integer codes.
  if (!is.object(data) && (tp == "double" || tp == "integer")) {
    as.vector(data)                 # strips names/attributes, keeps the type
  } else {
    as.vector(data, "numeric")      # logical, character, complex, list, factor
  }
}

#' Construct a ROIVolWindow without paying for the default initialize()
#'
#' Equivalent to \code{new("ROIVolWindow", data, space = space, ...)} --
#' \code{identical()} to it -- but roughly an order of magnitude cheaper,
#' which matters because searchlights build one of these per voxel.
#'
#' @keywords internal
#' @noRd
.new_roi_vol_window <- function(data, space, coords, center_index, parent_index) {
  data <- .as_roi_data(data)

  # The invariants the class validity function checks...
  if (!is.matrix(coords) || ncol(coords) != 3L) {
    stop("coords slot must be a matrix with 3 columns")
  }
  if (!is.vector(data)) {
    stop("'data' must be a vector")
  }
  if (length(data) != nrow(coords)) {
    stop("length of data vector must equal 'nrow(coords)'")
  }
  # ...plus the slot types that `@<-` would otherwise check on each assignment.
  if (!is(space, "NeuroSpace")) {
    stop("space slot must be a NeuroSpace object")
  }
  if (!is.integer(center_index) || !is.integer(parent_index)) {
    stop("center_index and parent_index must be integer")
  }

  # slot(check = FALSE) skips the per-assignment class lookup that `@<-` does.
  # Everything it would have checked is asserted above, once, and the result is
  # identical() to new() -- verified in test-roi-series-fastpaths.R.
  obj <- .roi_prototype("ROIVolWindow")
  slot(obj, ".Data", check = FALSE) <- data
  slot(obj, "space", check = FALSE) <- space
  slot(obj, "coords", check = FALSE) <- coords
  slot(obj, "center_index", check = FALSE) <- center_index
  slot(obj, "parent_index", check = FALSE) <- parent_index
  obj
}


# ---------------------------------------------------------------------------
# Shared semantics for the windowed ROI builders
#
# spherical_roi(), cuboid_roi() and square_roi() used to each implement centroid
# validation, `nonzero` filtering and the centre/parent bookkeeping separately,
# and disagreed with one another on all three. The two helpers below define it
# once so the three builders cannot drift apart again.
# ---------------------------------------------------------------------------

#' Validate and normalise an ROI centroid
#'
#' Accepts a length-3 vector or a 1x3 matrix, requires it to lie inside the
#' volume, and returns it as an integer voxel coordinate. Non-integer input is
#' truncated (as \code{as.integer()} does) so that the returned value is the one
#' actually used both to build the neighbourhood and to locate the centre within
#' it -- previously the grid was built from the truncated value while the centre
#' was matched against the original, so a fractional centroid never matched.
#'
#' @keywords internal
#' @noRd
.roi_centroid <- function(centroid, vdim, what = "centroid") {
  if (is.matrix(centroid)) {
    if (ncol(centroid) != 3 || nrow(centroid) != 1) {
      cli::cli_abort("{.arg {what}} matrix must have exactly 1 row and 3 columns.")
    }
    centroid <- drop(centroid)
  }

  if (length(centroid) != 3) {
    cli::cli_abort("{.arg {what}} must have length 3, not {length(centroid)}.")
  }
  if (!is.numeric(centroid) || anyNA(centroid) || any(!is.finite(centroid))) {
    cli::cli_abort("{.arg {what}} must contain finite numeric values.")
  }
  if (!all(centroid >= 1)) {
    cli::cli_abort("All {.arg {what}} values must be >= 1.")
  }

  centroid <- as.integer(centroid)
  if (centroid[1] > vdim[1] || centroid[2] > vdim[2] || centroid[3] > vdim[3]) {
    cli::cli_abort("{.arg {what}} ({.val {centroid}}) exceeds volume dimensions ({.val {vdim[1:3]}}).")
  }

  centroid
}

#' Assemble a ROIVolWindow from a candidate grid
#'
#' Applies the \code{nonzero} filter and fills in the centre/parent bookkeeping
#' with one definition shared by every windowed ROI builder:
#'
#' \itemize{
#'   \item \code{nonzero = TRUE} means what it says. The centre voxel is dropped
#'     along with any other zero-valued voxel; it is not forced back in.
#'   \item \code{center_index} is the row of the centre voxel within
#'     \code{coords}, or \code{NA_integer_} when the centre is not part of the
#'     ROI. It is never a fallback row, which would silently point at an
#'     unrelated voxel.
#'   \item \code{parent_index} is always the linear index of the requested
#'     centroid in the parent space, independent of any filtering -- that is
#'     what the slot documents.
#'   \item \code{coords} is an integer matrix without dimnames.
#' }
#'
#' @keywords internal
#' @noRd
.build_roi_window <- function(bvol, centroid, grid, vals, nonzero) {
  bspace <- space(bvol)
  bdim <- dim(bspace)
  parent_index <- as.integer((centroid[3] - 1L) * bdim[1] * bdim[2] +
                             (centroid[2] - 1L) * bdim[1] + centroid[1])

  if (nrow(grid) > 0L && isTRUE(nonzero)) {
    keep <- vals != 0
    grid <- grid[keep, , drop = FALSE]
    vals <- vals[keep]
  }

  # The compiled sphere path already yields undecorated integer coords; the
  # expand.grid-based builders do not, so normalise only when needed.
  if (storage.mode(grid) != "integer") storage.mode(grid) <- "integer"
  if (!is.null(dimnames(grid))) dimnames(grid) <- NULL

  if (nrow(grid) == 0L) {
    return(.new_roi_vol_window(numeric(0), bspace, matrix(integer(0), ncol = 3, nrow = 0),
                               NA_integer_, parent_index))
  }

  center_index <- which(grid[, 1] == centroid[1] &
                        grid[, 2] == centroid[2] &
                        grid[, 3] == centroid[3])
  center_index <- if (length(center_index) == 0L) NA_integer_ else as.integer(center_index[1])

  .new_roi_vol_window(vals, bspace, grid, center_index, parent_index)
}


#' ROI Coordinates
#'
#' @title Create ROI Coordinates Object
#' @description
#' Creates an \code{\linkS4class{ROICoords}} object from a matrix of coordinates
#' representing points in 3D space.
#'
#' @param coords A matrix with 3 columns representing (x, y, z) coordinates
#'
#' @return An \code{\linkS4class{ROICoords}} object
#'
#' @examples
#' coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
#' roi_coords <- ROICoords(coords)
#'
#' @export
ROICoords <- function(coords) {
  if (!is.matrix(coords) || ncol(coords) != 3) {
    cli::cli_abort("{.arg coords} must be a matrix with 3 columns (x,y,z).")
  }
  new("ROICoords", coords=coords)
}

#' ROI Volume
#'
#' @title Create ROI Volume Object
#' @description
#' Creates an \code{\linkS4class{ROIVol}} object representing a set of values
#' at specific 3D coordinates within a spatial reference system.
#'
#' @param space A \code{\linkS4class{NeuroSpace}} object defining the spatial reference
#' @param coords A matrix with 3 columns representing (x,y,z) coordinates
#' @param data A numeric vector of values corresponding to each coordinate
#'
#' @return An \code{\linkS4class{ROIVol}} object
#'
#' @examples
#' space <- NeuroSpace(c(64,64,64))
#' coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
#' data <- c(1.5, 2.5)
#' roi_vol <- ROIVol(space, coords, data)
#'
#' @export
ROIVol <- function(space, coords, data) {
  if (!is.matrix(coords) || ncol(coords) != 3) {
    cli::cli_abort("{.arg coords} must be a matrix with 3 columns (x,y,z).")
  }
  if (length(data) != nrow(coords)) {
    cli::cli_abort("Length of {.arg data} ({length(data)}) must match number of coordinates ({nrow(coords)}).")
  }
  new("ROIVol", data, space=space, coords=coords)
}


#' Extract or replace parts of an object
#' @name [
#' @rdname extract-methods
#' @aliases [,NeuroVol,ROICoords,missing,ANY-method
#'          [,ROICoords,numeric,missing,ANY-method
#'          [,ROIVol,ROICoords,missing,ANY-method
#'          [,ROIVol,ROICoords,numeric,ANY-method
#'          [,ROIVol,logical,missing,ANY-method
#'          [,ROIVol,logical,numeric,ANY-method
#'          [,ROIVol,matrix,missing,ANY-method
#'          [,ROIVol,matrix,numeric,ANY-method
#'          [,ROIVol,missing,missing,ANY-method
#'          [,ROIVol,missing,numeric,ANY-method
#'          [,ROIVol,numeric,numeric,ANY-method
#'          [,SparseNeuroVol,numeric,numeric,ANY-method
#'          [,ROIVol,numeric,missing,ANY-method
#'          [,NeuroVol,ROIVol,missing,ANY-method
#'          [,DenseNeuroVol,ROIVol,missing,ANY-method
#'          [,DenseNeuroVol,integer,missing,ANY-method
#'          [,DenseNeuroVol,numeric,missing,ANY-method
#' @param x The object to extract from
#' @param i Index specifying elements to extract
#' @param j Second index (if applicable)
#' @param k Third index for 3D objects (if applicable)
#' @param ... Additional arguments passed to methods
#' @param drop Whether to drop dimensions of length 1
#' @return A subset of the input object, with dimensions depending on the indexing and the `drop` parameter.
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "numeric", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            ROIVol(space(x), x@coords[i,,drop=FALSE], x@.Data[i])
          }
)


#' @export
#' @rdname extract-methods
setMethod(f="[", signature=signature(x = "ROIVol", i = "logical", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            ROIVol(space(x), x@coords[i,,drop=FALSE], x@.Data[i])
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "ROICoords", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            callGeneric(x, i@coords)
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "matrix", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            stopifnot(ncol(i) == 3)
            stopifnot(nrow(i) == nrow(x@coords))

            ROIVol(space(x), i, x@.Data)
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "missing", j = "missing"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            x@.Data
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "missing", j = "numeric"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            x@coords[,j,drop=FALSE]
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "numeric", j = "numeric"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            x@coords[i,j,drop=FALSE]
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "logical", j = "numeric"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            x@coords[i,j,drop=FALSE]
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "ROICoords", j = "numeric"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            callGeneric(x, i@coords, j)
          }
)


#' @rdname extract-methods
#' @export
setMethod(f="[", signature=signature(x = "ROIVol", i = "matrix", j = "numeric"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            stopifnot(ncol(i) == 3)
            stopifnot(nrow(i) == nrow(x@coords))

            i[,j,drop=FALSE]
          }
)


#' @export
#' @rdname dim-methods
setMethod(f="dim", signature=signature(x = "ROIVol"),
          def=function(x) {
            c(nrow(x@coords), ncol(x@coords))
          })



#' @export
#' @rdname length-methods
setMethod(f="length", signature=signature(x = "ROIVol"),
          def=function(x) {
            nrow(x@coords)
          })


#' @export
#' @rdname indices-methods
setMethod(f="indices", signature=signature(x = "ROIVol"),
          def=function(x) {
            grid_to_index(x@space, x@coords)
          })



#' @rdname as.numeric-methods
#' @export
#' @return A numeric vector of length \code{nrow(x@coords)}
setMethod(f="as.numeric", signature=signature(x = "ROIVol"),
          def=function(x) {
            x@.Data
          })


#' @rdname as.sparse-methods
#' @export
setMethod(f="as.sparse", signature=signature(x = "ROIVol"),
          def=function(x) {
            SparseNeuroVol(x@.Data, space(x), indices=indices(x))
          })


#' Convert object to logical type
#' @name as.logical
#' @rdname as.logical-methods
#' @aliases as.logical,ROIVol-method
#' @export
setMethod(f="as.logical", signature=signature(x = "ROIVol"),
          def=function(x) {
            ret <- array(FALSE, dim(x@space))
            ret[indices(x)] <- TRUE
            LogicalNeuroVol(ret, space(x))
          })



#' @rdname coords-methods
#' @param real Logical indicating whether to return real coordinates
#' @return A matrix of coordinates
#' @export
setMethod(f="coords", signature=signature(x = "ROIVol"),
          def=function(x, real=FALSE) {
            if(real) {
              input <- t(cbind(x@coords - 0.5, rep(1, nrow(x@coords))) );
              ret <- t(trans(x@space) %*% input);
              ret[,1:3, drop=FALSE]
            } else {
              x@coords
            }
          })


#' @export
#' @rdname show-methods
setMethod("show", "ROIVol", function(object) {
  show_header("ROIVol", paste(length(object@.Data), "voxels"))
  show_field("Coords", paste(nrow(coords(object)), "x 3"))
  rng <- range(object@.Data, na.rm = TRUE)
  show_field("Range", sprintf("[%.3f, %.3f]", rng[1], rng[2]))
})


#' @rdname show-methods
#' @export
setMethod(f="show", signature=signature(object = "ROICoords"),
          def=function(object) {
            cat("ROICoords\n")
            cat("  Dimension      :", dim(object@coords), "\n")
          })


#' @export
#' @rdname dim-methods
#' @return A numeric vector of length 2 containing the dimensions of the ROICoords object.
setMethod(f="dim", signature=signature(x = "ROICoords"),
          def=function(x) {
            dim(x@coords)
          })


#' @export
#' @rdname length-methods
#' @return An integer representing the number of coordinates in the ROICoords object.
setMethod(f="length", signature=signature(x = "ROICoords"),
          def=function(x) {
            nrow(x@coords)
          })


#' @rdname extract-methods
#' @export
setMethod("[", signature(x="ROICoords", i="numeric", j="missing", drop="ANY"),
          function(x, i, j, ..., drop) {
            new("ROICoords", coords=x@coords[i,,drop=FALSE])
          })



#' @rdname extract-methods
#' @export
setMethod("[", signature(x="ROIVol", i="numeric", j="missing", drop="ANY"),
          function (x, i, j, ..., drop) {
            ROIVol(x@space, x@coords[i,,drop=FALSE], x@.Data[i])
          })

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
ROIVec <- function(vspace, coords, data = matrix(1, nrow = 1, ncol = nrow(coords))) {
  # Accept vector input and reshape to 1 x ncoords; otherwise ensure columns match coords
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  if (ncol(data) != nrow(coords)) {
    stop("ROIVec: number of columns in 'data' must equal number of rows in 'coords'")
  }
  new("ROIVec", space = vspace, coords = coords, data)
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

  # expand.grid() varies its FIRST argument fastest, so listing z first and then
  # reordering the columns yields rows in (x, y, z) lexicographic order -- the
  # same order the compiled spherical path emits. Every ROI builder therefore
  # returns coordinates in one canonical order, at no cost.
  if (fixdim == 3) {
    grid <- as.matrix(expand.grid(z=centroid[3], y=coords[[2]], x=coords[[1]]))
  } else if (fixdim == 2) {
    grid <- as.matrix(expand.grid(z=coords[[2]], y=centroid[2], x=coords[[1]]))
  } else if (fixdim == 1) {
    grid <- as.matrix(expand.grid(z=coords[[2]], y=coords[[1]], x=centroid[1]))
  }

  grid[, c("x", "y", "z"), drop = FALSE]

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

  # See .makeSquareGrid(): reversed argument order gives (x, y, z) lexicographic
  # rows, matching every other ROI builder.
  grid <- as.matrix(expand.grid(z=coords[[3]], y=coords[[2]], x=coords[[1]]))
  grid[, c("x", "y", "z"), drop = FALSE]
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
  centroid <- .roi_centroid(centroid, dim(bvol))

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

  .build_roi_window(bvol, centroid, grid, vals, nonzero)
}


#' Create A Cuboid Region of Interest
#'
#' @param bvol an \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroid the center of the cube in \emph{voxel} coordinates
#' @param surround the number of voxels on either side of the central voxel. A \code{vector} of length 3.
#' @param fill optional value(s) to assign to data slot.
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{NeuroSpace} then this argument is ignored.
#' @return An instance of class \code{ROIVol} representing the cuboid region of interest, containing the coordinates and values of voxels within the specified region.
#' @examples
#'  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
#'  cube <- cuboid_roi(sp1, c(5,5,5), 3)
#'  vox <- coords(cube)
#'  cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=5)
#'
#'
#' @export
cuboid_roi <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE) {
  centroid <- .roi_centroid(centroid, dim(bvol))

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

  .build_roi_window(bvol, centroid, grid, vals, nonzero)

}

#' Cached spherical offset template
#'
#' The voxel offsets forming a spherical neighbourhood depend only on the
#' radius and the voxel spacing, so they are computed once per
#' \code{(radius, spacing)} pair and reused across centres. This is what makes
#' an exhaustive searchlight cheap: per centre the work drops to translating
#' the template and clipping it to the volume bounds.
#'
#' @keywords internal
#' @noRd
.sphere_offset_cache <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.sphere_offsets <- function(radius, spacing) {
  if (length(radius) != 1L) {
    stop("'radius' must be a single value")
  }
  spacing <- as.numeric(spacing)[1:3]
  key <- paste0(sprintf("%.17g", radius), "|",
                paste(sprintf("%.17g", spacing), collapse = ","))

  hit <- .sphere_offset_cache[[key]]
  if (!is.null(hit)) {
    return(hit)
  }

  # Bound the cache: analyses sweep a handful of radii, not thousands.
  if (length(ls(.sphere_offset_cache, all.names = TRUE)) >= 64L) {
    rm(list = ls(.sphere_offset_cache, all.names = TRUE), envir = .sphere_offset_cache)
  }

  off <- sphere_offsets_cpp(radius, spacing)
  assign(key, off, envir = .sphere_offset_cache)
  off
}

#' Voxel coordinates of a spherical neighbourhood, 1-based integers
#'
#' @param use_cpp Deprecated and ignored. The former \code{FALSE} branch used
#'   \code{dbscan::frNN} over a candidate cube sized with \code{round()} rather
#'   than \code{ceil()}, and compared distances exclusively at the boundary. It
#'   therefore returned neighbourhoods that were both offset by one voxel and
#'   missing boundary voxels: checked against brute-force enumeration over 384
#'   geometries, the compiled path was correct in all of them and the fallback
#'   wrong in 35. Passing \code{FALSE} now warns and uses the compiled path.
#'
#' @keywords internal
#' @noRd
make_spherical_grid <- function(bvol, centroid, radius, use_cpp = TRUE) {

  if (!isTRUE(use_cpp)) {
    warning("'use_cpp = FALSE' is deprecated and ignored: the dbscan fallback ",
            "returned incorrect neighbourhoods. The compiled path is always used.",
            call. = FALSE)
  }

  vspacing <- spacing(bvol)

  if (radius < min(vspacing)) {
    stop("'radius' is too small; must be greater than at least one voxel dimension in image")
  }

  vdim <- dim(bvol)
  centroid <- as.integer(centroid)

  # Translate the cached offset template to this centre and clip to bounds.
  # Equivalent to local_sphere() -- same membership test, same emission order --
  # but the sphere itself is built once per (radius, spacing). Emitted 1-based
  # so callers need no arithmetic on the result; the 0-based convention only
  # existed for the retired dbscan branch.
  sphere_at_cpp(.sphere_offsets(radius, vspacing),
                as.integer(centroid[1:3]),
                as.integer(vdim[1:3]),
                base0 = FALSE)
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
#' @param use_cpp Deprecated and ignored; the compiled implementation is always
#'   used. Passing \code{FALSE} warns. The former fallback returned
#'   neighbourhoods offset by one voxel and missing boundary voxels.
#' @return an instance of class \code{ROIVol}
#' @seealso [spherical_roi_set()] for efficiently creating many spherical ROIs,
#'   [series_roi()] and [coords()] for extracting time series and coordinates from ROIs,
#'   and the vignette: \code{vignette("regionOfInterest", package = "neuroim2")}.
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
#'  ## Create multiple spherical ROIs at once (preferred):
#'  centers <- rbind(c(5,5,5), c(3,3,3), c(7,7,7))
#'  vols <- spherical_roi_set(bvol = sp1,
#'                           centroids = centers, radius = 3.5, fill = 1)
#'  length(vols)  # 3
#'
#'  ## Equivalent, less efficient lapply variant:
#'  vols2 <- lapply(seq_len(nrow(centers)), function(i) {
#'    spherical_roi(sp1, centers[i,], radius = 3.5, fill = 1)
#'  })
#'
#'  # create an ROI centered around the real-valued coordinates: x=5, y=5, z=5
#'  vox <- coord_to_grid(sp1, c(5, 5, 5))
#'  cube <- spherical_roi(sp1, vox, 3.5)
#' @export
spherical_roi <- function(bvol, centroid, radius, fill=NULL, nonzero=FALSE, use_cpp=TRUE) {
  centroid <- .roi_centroid(centroid, dim(bvol))

  if (is.null(fill) && is(bvol, "NeuroSpace")) {
    fill <- 1
  }

  # 1-based integer coords, already ordered by (x, y, z).
  grid <- make_spherical_grid(bvol, centroid, radius, use_cpp=use_cpp)

  if (nrow(grid) == 0) {
    return(.build_roi_window(bvol, centroid, grid, numeric(0), nonzero))
  }

  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    as.numeric(bvol[grid])
  }

  .build_roi_window(bvol, centroid, grid, vals, nonzero)
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



# Coerce ROIVec to matrix.
setAs(from="ROIVec", to="matrix", function(from) {
  ind <- indices(from)
  roi_vector_matrix(from@.Data, refspace=from@space, indices=ind,
                    coords=index_to_coord(drop_dim(from@space),
                                          as.numeric(ind)))

})

# Coerce ROIVol to DenseNeuroVol.
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
            cds = coords(x)
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
#' @aliases values,ROIVec-method
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
#' @rdname show-methods
setMethod("show", "ROIVec", function(object) {
  np <- ncol(object@.Data)
  nf <- nrow(object@.Data)
  show_header("ROIVec", paste(np, "voxels x", nf, "features"))
  show_rule("Structure")
  show_field("Voxels", np)
  show_field("Features", nf)
  show_field("Size", format_mem(object))
  show_rule("Data")
  rng <- range(object@.Data, na.rm = TRUE)
  show_field("Range", sprintf("[%.3f, %.3f]", rng[1], rng[2]))
  nas <- sum(is.na(object@.Data))
  if (nas > 0) show_field("NAs", format(nas, big.mark = ","))
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
  cvoxel.ind <- t(t(voxel.ind) * vdim)

  ## the coordinates ofthe voxels (i.e. after multiplying by pixel dims)
  coords <- t(apply(cvoxel.ind, 1, function(vals) sign(vals)* ifelse(vals == 0, 0, abs(vals)-.5)))

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



#' Create Multiple Spherical Regions of Interest
#'
#' @description
#' This function generates multiple spherical ROIs simultaneously, centered at the provided
#' voxel coordinates. It is more efficient than calling \code{spherical_roi} multiple times
#' when you need to create many ROIs.
#'
#' @param bvol A \code{NeuroVol} or \code{NeuroSpace} instance
#' @param centroids A matrix of voxel coordinates where each row represents a centroid (i,j,k)
#' @param radius The radius in real units (e.g. millimeters) of the spherical ROIs
#' @param fill Optional value(s) to store as data. If provided, must be either a single value
#'   or a vector with length equal to the number of ROIs
#' @param nonzero If \code{TRUE}, keep only nonzero elements from \code{bvol}
#'
#' @return A list of \code{ROIVolWindow} objects, one for each centroid
#'
#' @examples
#' # Create a NeuroSpace object
#' sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))
#'
#' # Create multiple ROIs centered at different voxel coordinates
#' centroids <- matrix(c(5,5,5, 3,3,3, 7,7,7), ncol=3, byrow=TRUE)
#' rois <- spherical_roi_set(sp1, centroids, 3.5)
#'
#' # Create ROIs with specific fill values
#' rois <- spherical_roi_set(sp1, centroids, 3.5, fill=c(1,2,3))
#'
#' @export
spherical_roi_set <- function(bvol, centroids, radius, fill=NULL, nonzero=FALSE) {
  # Pre-allocate result list
  n <- nrow(centroids)
  result_list <- vector("list", n)

  # Just loop and call spherical_roi() for each centroid
  for (i in seq_len(n)) {
    result_list[[i]] <- spherical_roi(
      bvol = bvol,
      centroid = centroids[i,],
      radius = radius,
      fill = if (!is.null(fill)) { if (length(fill) == 1) fill else fill[i] } else NULL,
      nonzero = nonzero,
      use_cpp = TRUE
    )
  }

  result_list
}

# spherical_roi_set <- function(bvol, centroids, radius, fill=NULL, nonzero=FALSE) {
#   if (!is.matrix(centroids)) {
#     stop("centroids must be a matrix with 3 columns")
#   }
#   if (ncol(centroids) != 3) {
#     stop("centroids must have exactly 3 columns (i,j,k coordinates)")
#   }
#
#   bspace <- space(bvol)
#   vspacing <- spacing(bvol)
#   vdim <- dim(bvol)
#
#   if (is.null(fill) && is(bvol, "NeuroSpace")) {
#     fill <- 1
#   }
#
#   if (!is.null(fill)) {
#     if (length(fill) != 1 && length(fill) != nrow(centroids)) {
#       stop("fill must be either length 1 or match the number of ROIs")
#     }
#   }
#
#   # Call local_spheres returning zero-based coords
#   sphere_coords_list <- local_spheres(
#     centers = centroids,
#     radius = radius,
#     spacing = vspacing,
#     dim = vdim
#   )
#
#   result_list <- vector("list", nrow(centroids))
#
#   for (i in seq_len(nrow(centroids))) {
#     sphere_coords <- sphere_coords_list[[i]]
#
#     # Convert to integer and now adjust +1 to move from zero-based to one-based
#     sphere_coords <- as.matrix(sphere_coords)
#     mode(sphere_coords) <- "integer"
#     sphere_coords <- sphere_coords + 1L
#
#     sphere_idx <- grid_to_index(bspace, sphere_coords)
#
#     if (nonzero) {
#       keep <- bvol[sphere_idx] != 0
#       sphere_coords <- sphere_coords[keep, , drop=FALSE]
#       sphere_idx <- sphere_idx[keep]
#     }
#
#     roi_fill <- if (!is.null(fill)) {
#       if (length(fill) == 1) fill else fill[i]
#     } else {
#       as.numeric(bvol[sphere_coords])
#     }
#
#     # Center index calculation remains the same
#     center_coords <- matrix(centroids[i,], ncol=3)
#     mode(center_coords) <- "integer"
#     parent_idx <- grid_to_index(bspace, center_coords)[1]
#
#     center_index <- which(
#       sphere_coords[,1] == centroids[i,1] &
#         sphere_coords[,2] == centroids[i,2] &
#         sphere_coords[,3] == centroids[i,3]
#     )[1]
#
#     if (is.na(center_index)) center_index <- 1L
#
#     result_list[[i]] <- new("ROIVolWindow",
#                             roi_fill,
#                             space=bspace,
#                             coords=sphere_coords,
#                             center_index=as.integer(center_index),
#                             parent_index=as.integer(parent_idx))
#   }
#
#   result_list
# }


#' @rdname extract-methods
#' @export
setMethod("[", signature(x="NeuroVol", i="ROICoords", j="missing", drop="ANY"),
          function(x, i, j, ..., drop) {
            callGeneric(x, i@coords, drop=drop)
          })

 
