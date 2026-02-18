#' @include all_generic.R
NULL
#' @include all_class.R
NULL

#' Pre-defined null axis
#' @export
None <- new("NamedAxis", axis="None")

#' Pre-defined null axis set
#' @export
NullAxis <- new("AxisSet", ndim=as.integer(0))

#' Pre-defined anatomical axes
#'
#' These constants define standard anatomical axes used in neuroimaging.
#' Each axis has a defined direction vector in 3D space.
#'
#' @name anatomical_axes
#' @export
LEFT_RIGHT <- new("NamedAxis", axis="Left-to-Right", direction=c(1,0,0))

#' @rdname anatomical_axes
#' @export
RIGHT_LEFT <- new("NamedAxis", axis="Right-to-Left", direction=c(-1,0,0))

#' @rdname anatomical_axes
#' @export
ANT_POST   <- new("NamedAxis", axis="Anterior-to-Posterior", direction=c(0,-1,0))

#' @rdname anatomical_axes
#' @export
POST_ANT   <- new("NamedAxis", axis="Posterior-to-Anterior", direction=c(0,1,0))

#' @rdname anatomical_axes
#' @export
INF_SUP    <- new("NamedAxis", axis="Inferior-to-Superior", direction=c(0,0,1))

#' @rdname anatomical_axes
#' @export
SUP_INF    <- new("NamedAxis", axis="Superior-to-Inferior", direction=c(0,0,-1))

#' Match anatomical axis abbreviation to full axis object
#'
#' @param firstAxis Character string specifying the axis. Can be full name
#'   (e.g. "LEFT") or abbreviation (e.g. "L")
#' @return A NamedAxis object corresponding to the specified axis
#' @keywords internal
#' @noRd
matchAxis <- function(firstAxis) {
  switch(toupper(firstAxis),
         "LEFT" = LEFT_RIGHT,
         "L" = LEFT_RIGHT,
         "RIGHT" = RIGHT_LEFT,
         "R" = RIGHT_LEFT,
         "ANTERIOR" = ANT_POST,
         "A" = ANT_POST,
         "POSTERIOR" = POST_ANT,
         "P" = POST_ANT,
         "INFERIOR" = INF_SUP,
         "I" = INF_SUP,
         "SUPERIOR" = SUP_INF,
         "S" = SUP_INF)
}

#' Time axis
#' @description Represents the temporal dimension in neuroimaging data
#' @export
TIME <- new("NamedAxis", axis="Time")

#' Time axis set
#' @description A one-dimensional axis set representing time
#' @export
TimeAxis <- new("AxisSet1D", ndim=as.integer(1), i=TIME)

#' Create a one-dimensional axis set
#'
#' @param i A NamedAxis object representing the axis
#' @return An AxisSet1D object
#' @keywords internal
#' @noRd
AxisSet1D <- function(i) {
  new("AxisSet1D", ndim = as.integer(1), i = i)
}

#' Create a two-dimensional axis set
#'
#' @param i A NamedAxis object representing the first axis
#' @param j A NamedAxis object representing the second axis
#' @return An AxisSet2D object
#' @keywords internal
#' @noRd
AxisSet2D <- function(i, j) {
	new("AxisSet2D", ndim = as.integer(2), i = i, j = j)
}

#' Create a three-dimensional axis set
#'
#' @param i A NamedAxis object representing the first axis
#' @param j A NamedAxis object representing the second axis
#' @param k A NamedAxis object representing the third axis
#' @return An AxisSet3D object
#' @keywords internal
#' @noRd
AxisSet3D <- function(i, j, k) {
	new("AxisSet3D", ndim = as.integer(3), i = i, j = j, k = k)
}

#' Get permutation matrix from axis set
#'
#' @param x An AxisSet2D object
#' @param ... Additional arguments (not used)
#' @return A matrix representing the axis directions
#' @export
setMethod(f="perm_mat", signature=signature(x = "AxisSet2D"),
          def=function(x, ...) {
            cbind(x@i@direction, x@j@direction)
          })

#' Get permutation matrix from axis set
#'
#' @param x An AxisSet3D object
#' @param ... Additional arguments (not used)
#' @return A matrix representing the axis directions
#' @export
setMethod(f="perm_mat", signature=signature(x = "AxisSet3D"),
          def=function(x, ...) {
            cbind(x@i@direction, x@j@direction, x@k@direction)
          })

#' Get permutation matrix from NeuroSpace
#'
#' @param x A NeuroSpace object
#' @param ... Additional arguments (not used)
#' @return A numeric N x N \code{matrix} representing the permutation transform, where N is the dimensionality of the image.
#' @export
#' @rdname perm_mat-methods
setMethod(f="perm_mat", signature=signature(x = "NeuroSpace"),
          def=function(x, ...) {
            callGeneric(x@axes)
          })


#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x = "AxisSet2D", dimnum="numeric"),
          def=function(x, dimnum) {
            stopifnot(length(dimnum) == 1)
            if (dimnum == 1) {
              AxisSet1D(x@j)
            } else if (dimnum == 2) {
              AxisSet1D(x@i)
            } else {
              stop(paste("illegal dimnum: ",dimnum, "for axis with 2 dimensions"))
            }
          })


#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x = "AxisSet2D", dimnum="missing"),
          def=function(x, dimnum) {
            AxisSet1D(x@i)
          })

#' @rdname drop_dim-methods
#' @export
setMethod(f="drop_dim", signature=signature(x = "AxisSet3D", dimnum="numeric"),
          def=function(x, dimnum) {
            stopifnot(length(dimnum) == 1)
            if (dimnum == 1) {
              AxisSet2D(x@j, x@k)
            } else if (dimnum == 2) {
              AxisSet2D(x@i, x@k)
            } else if (dimnum == 3) {
              AxisSet2D(x@i, x@j)
            } else {
              stop(paste("illegal dimnum: ", dimnum, " for axis with 2 dimensions"))
            }
          })

#' Drop dimension from axis set
#'
#' @param x An AxisSet3D object
#' @param dimnum Numeric index of dimension to drop (optional)
#' @export
#' @rdname drop_dim-methods
setMethod(f="drop_dim", signature=signature(x = "AxisSet3D", dimnum="missing"),
           def=function(x, dimnum) {
             AxisSet2D(x@i, x@j)
           })

#' Get number of dimensions in axis set
#'
#' @param x An AxisSet object
#' @param ... Additional arguments (not used)
#' @return An integer representing the number of dimensions in \code{x}.
#' @export
setMethod(f="ndim",signature(x= "AxisSet"), function(x, ...) { x@ndim })

#' Print method for NamedAxis objects
#'
#' @param x A NamedAxis object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis
#' @keywords internal
#' @noRd
setMethod(f="print_", signature=signature("NamedAxis"),
    def=function(x, ...) {
        x@axis
    })

#' Show method for NamedAxis objects
#'
#' @param object A NamedAxis object
#' @return Invisibly returns \code{NULL}, called for its side effect of displaying the object.
#' @export
#' @rdname show-methods
setMethod(f="show", signature("NamedAxis"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("NamedAxis"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("-", 30), collapse="")), "\n")
        cat(crayon::white(object@axis), "\n")
    })

#' @export
#' @rdname show-methods
setMethod(f="show", signature=signature("AxisSet1D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet1D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("-", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis:"), crayon::white(object@i@axis), "\n")
    })

#' Print method for AxisSet2D objects
#'
#' @param x An AxisSet2D object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis set
#' @keywords internal
#' @noRd
setMethod(f="print_", signature=signature("AxisSet2D"),
    def=function(x, ...) {
        paste(x@i@axis, "x", x@j@axis)
    })

#' @rdname show-methods
#' @export
setMethod(f="show", signature=signature("AxisSet2D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet2D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("-", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis 1:"), crayon::white(object@i@axis), "\n")
        cat(crayon::yellow("Axis 2:"), crayon::white(object@j@axis), "\n")
    })

#' Print method for AxisSet3D objects
#'
#' @param x An AxisSet3D object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis set
#' @keywords internal
#' @noRd
setMethod(f="print_", signature=signature("AxisSet3D"),
    def=function(x, ...) {
        paste(x@i@axis, "x", x@j@axis, "x", x@k@axis)
    })

#' @rdname show-methods
#' @export
setMethod(f="show", signature=signature("AxisSet3D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet3D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("-", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis 1:"), crayon::white(object@i@axis), "\n")
        cat(crayon::yellow("Axis 2:"), crayon::white(object@j@axis), "\n")
        cat(crayon::yellow("Axis 3:"), crayon::white(object@k@axis), "\n")
    })

#' @rdname show-methods
#' @export
setMethod(f="show", signature=signature("AxisSet4D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet4D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("-", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis 1:"), crayon::white(object@i@axis), "\n")
        cat(crayon::yellow("Axis 2:"), crayon::white(object@j@axis), "\n")
        cat(crayon::yellow("Axis 3:"), crayon::white(object@k@axis), "\n")
        cat(crayon::yellow("Axis 4:"), crayon::white(object@t@axis), "\n")
    })

#' Pre-defined 2D orientation configurations
#'
#' A list of standard 2D anatomical orientations used in neuroimaging.
#' Each orientation defines a pair of anatomical axes.
#'
#' @export
OrientationList2D <- list(
	SAGITTAL_AI = AxisSet2D(ANT_POST, INF_SUP),
	SAGITTAL_PI = AxisSet2D(POST_ANT, INF_SUP),
	SAGITTAL_PS = AxisSet2D(POST_ANT, SUP_INF),
	SAGITTAL_AS = AxisSet2D(ANT_POST, SUP_INF),
	SAGITTAL_IA = AxisSet2D(INF_SUP, ANT_POST),
	SAGITTAL_IP = AxisSet2D(INF_SUP, POST_ANT),
	SAGITTAL_SP = AxisSet2D(SUP_INF, POST_ANT),
	SAGITTAL_SA = AxisSet2D(SUP_INF, ANT_POST),

	CORONAL_LI = AxisSet2D(LEFT_RIGHT, INF_SUP),
	CORONAL_RI = AxisSet2D(RIGHT_LEFT, INF_SUP),
	CORONAL_RS = AxisSet2D(RIGHT_LEFT, SUP_INF),
	CORONAL_LS = AxisSet2D(LEFT_RIGHT, SUP_INF),
	CORONAL_IL = AxisSet2D(INF_SUP, LEFT_RIGHT),
	CORONAL_IR = AxisSet2D(INF_SUP, RIGHT_LEFT),
	CORONAL_SR = AxisSet2D(SUP_INF, RIGHT_LEFT),
	CORONAL_SL = AxisSet2D(SUP_INF, LEFT_RIGHT),


	AXIAL_LA = AxisSet2D(LEFT_RIGHT, ANT_POST),
	AXIAL_RA = AxisSet2D(RIGHT_LEFT, ANT_POST),
	AXIAL_RP = AxisSet2D(RIGHT_LEFT, POST_ANT),
	AXIAL_LP = AxisSet2D(LEFT_RIGHT, POST_ANT),
	AXIAL_AL = AxisSet2D(ANT_POST, LEFT_RIGHT),
	AXIAL_AR = AxisSet2D(ANT_POST, RIGHT_LEFT),
	AXIAL_PL = AxisSet2D(POST_ANT, LEFT_RIGHT),
	AXIAL_PR = AxisSet2D(POST_ANT, RIGHT_LEFT))

#' Pre-defined 3D orientation configurations
#'
#' A list of standard 3D anatomical orientations used in neuroimaging.
#' Each orientation defines a triplet of anatomical axes.
#'
#' @export
OrientationList3D <- list(
	SAGITTAL_AIL = AxisSet3D(ANT_POST, INF_SUP, LEFT_RIGHT),
	SAGITTAL_PIL = AxisSet3D(POST_ANT, INF_SUP, LEFT_RIGHT),
	SAGITTAL_PSL = AxisSet3D(POST_ANT, SUP_INF, LEFT_RIGHT),
	SAGITTAL_ASL = AxisSet3D(ANT_POST, SUP_INF, LEFT_RIGHT),
	SAGITTAL_IAL = AxisSet3D(INF_SUP,  ANT_POST, LEFT_RIGHT),
	SAGITTAL_IPL = AxisSet3D(INF_SUP,  POST_ANT, LEFT_RIGHT),
	SAGITTAL_SPL = AxisSet3D(SUP_INF,  POST_ANT, LEFT_RIGHT),
	SAGITTAL_SAL = AxisSet3D(SUP_INF,  ANT_POST, LEFT_RIGHT),

	SAGITTAL_AIR = AxisSet3D(ANT_POST, INF_SUP, RIGHT_LEFT),
	SAGITTAL_PIR = AxisSet3D(POST_ANT, INF_SUP, RIGHT_LEFT),
	SAGITTAL_PSR = AxisSet3D(POST_ANT, SUP_INF, RIGHT_LEFT),
	SAGITTAL_ASR = AxisSet3D(ANT_POST, SUP_INF, RIGHT_LEFT),
	SAGITTAL_IAR = AxisSet3D(INF_SUP,  ANT_POST, RIGHT_LEFT),
	SAGITTAL_IPR = AxisSet3D(INF_SUP,  POST_ANT, RIGHT_LEFT),
	SAGITTAL_SPR = AxisSet3D(SUP_INF,  POST_ANT, RIGHT_LEFT),
	SAGITTAL_SAR = AxisSet3D(SUP_INF,  ANT_POST, RIGHT_LEFT),

	CORONAL_LIA = AxisSet3D(LEFT_RIGHT, INF_SUP, ANT_POST),
	CORONAL_RIA = AxisSet3D(RIGHT_LEFT, INF_SUP, ANT_POST),
	CORONAL_RSA = AxisSet3D(RIGHT_LEFT, SUP_INF, ANT_POST),
	CORONAL_LSA = AxisSet3D(LEFT_RIGHT, SUP_INF, ANT_POST),
	CORONAL_ILA = AxisSet3D(INF_SUP,    LEFT_RIGHT, ANT_POST),
	CORONAL_IRA = AxisSet3D(INF_SUP,    RIGHT_LEFT, ANT_POST),
	CORONAL_SRA = AxisSet3D(SUP_INF,    RIGHT_LEFT, ANT_POST),
	CORONAL_SLA = AxisSet3D(SUP_INF,    LEFT_RIGHT, ANT_POST),

	CORONAL_LIP = AxisSet3D(LEFT_RIGHT, INF_SUP, POST_ANT),
	CORONAL_RIP = AxisSet3D(RIGHT_LEFT, INF_SUP, POST_ANT),
	CORONAL_RSP = AxisSet3D(RIGHT_LEFT, SUP_INF, POST_ANT),
	CORONAL_LSP = AxisSet3D(LEFT_RIGHT, SUP_INF, POST_ANT),
	CORONAL_ILP = AxisSet3D(INF_SUP,    LEFT_RIGHT, POST_ANT),
	CORONAL_IRP = AxisSet3D(INF_SUP,    RIGHT_LEFT, POST_ANT),
	CORONAL_SRP = AxisSet3D(SUP_INF,    RIGHT_LEFT, POST_ANT),
	CORONAL_SLP = AxisSet3D(SUP_INF,    LEFT_RIGHT, POST_ANT),


	AXIAL_LAI = AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP),
	AXIAL_RAI = AxisSet3D(RIGHT_LEFT, ANT_POST, INF_SUP),
	AXIAL_RPI = AxisSet3D(RIGHT_LEFT, POST_ANT, INF_SUP),
	AXIAL_LPI = AxisSet3D(LEFT_RIGHT, POST_ANT, INF_SUP),
	AXIAL_ALI = AxisSet3D(ANT_POST,   LEFT_RIGHT, INF_SUP),
	AXIAL_ARI = AxisSet3D(ANT_POST,   RIGHT_LEFT, INF_SUP),
	AXIAL_PLI = AxisSet3D(POST_ANT,   LEFT_RIGHT, INF_SUP),
	AXIAL_PRI = AxisSet3D(POST_ANT,   RIGHT_LEFT, INF_SUP),

	AXIAL_LAS = AxisSet3D(LEFT_RIGHT, ANT_POST, SUP_INF),
	AXIAL_RAS = AxisSet3D(RIGHT_LEFT, ANT_POST, SUP_INF),
	AXIAL_RPS = AxisSet3D(RIGHT_LEFT, POST_ANT, SUP_INF),
	AXIAL_LPS = AxisSet3D(LEFT_RIGHT, POST_ANT, SUP_INF),
	AXIAL_ALS = AxisSet3D(ANT_POST,   LEFT_RIGHT, SUP_INF),
	AXIAL_ARS = AxisSet3D(ANT_POST,   RIGHT_LEFT, SUP_INF),
	AXIAL_PLS = AxisSet3D(POST_ANT,   LEFT_RIGHT, SUP_INF),
	AXIAL_PRS = AxisSet3D(POST_ANT,   RIGHT_LEFT, SUP_INF))


#' Find 3D anatomical orientation from axis abbreviations
#'
#' Creates a 3D anatomical orientation from axis abbreviations.
#'
#' @param axis1 Character string for first axis (default: "L" for Left)
#' @param axis2 Character string for second axis (default: "P" for Posterior)
#' @param axis3 Character string for third axis (default: "I" for Inferior)
#' @return An AxisSet3D object representing the anatomical orientation
#' @export
#' @examples
#' # Create orientation with default LPI axes
#' orient <- findAnatomy3D()
#' # Create orientation with custom axes
#' orient <- findAnatomy3D("R", "A", "S")
findAnatomy3D <- function(axis1="L", axis2="P", axis3="I") {
  res <- lapply(list(axis1, axis2, axis3), function(x) {
    matchAxis(x)
  })

  matchAnatomy3D(res[[1]], res[[2]], res[[3]])

}

#' Find matching 2D anatomical orientation
#'
#' Searches through pre-defined 2D orientations to find a match for given axes.
#'
#' @param axis1 First NamedAxis object
#' @param axis2 Second NamedAxis object
#' @return Matching AxisSet2D object or NULL if no match found
#' @keywords internal
#' @noRd
matchAnatomy2D <- function(axis1, axis2) {
	for (orient in OrientationList2D) {
		if (identical(orient@i,axis1) && identical(orient@j,axis2)) {
			return(orient)
		}
	}

	stop("did not find legal matching anatomical orientation for axes: ",
			axis1, axis2)

}

#' Find matching 3D anatomical orientation
#'
#' Searches through pre-defined 3D orientations to find a match for given axes.
#'
#' @param axis1 First NamedAxis object
#' @param axis2 Second NamedAxis object
#' @param axis3 Third NamedAxis object
#' @return Matching AxisSet3D object or NULL if no match found
#' @keywords internal
#' @noRd
matchAnatomy3D <- function(axis1, axis2, axis3) {
	for (orient in OrientationList3D) {
		if (identical(orient@i,axis1) && identical(orient@j,axis2) && identical(orient@k, axis3)) {
			return(orient)
		}
	}

	stop("did not find legal matching anatomical orientation for axes: ",
			axis1, axis2, axis3)

}

#' Orientation labels used for axis-code conversion
#'
#' @return List of `(negative_end, positive_end)` label pairs.
#' @keywords internal
#' @noRd
.default_axcode_labels <- function() {
  list(c("L", "R"), c("P", "A"), c("I", "S"))
}

#' Validate and normalize orientation matrix
#'
#' Orientation matrices are `n x 2` where column 1 gives target axis index
#' (1-based, or `NA` for dropped axes) and column 2 gives direction (`-1` or `1`).
#'
#' @param ornt Orientation matrix-like object.
#' @param allow_dropped Whether `NA` axes are allowed.
#' @param arg_name Name for error messages.
#' @return Numeric matrix with columns `axis` and `flip`.
#' @keywords internal
#' @noRd
.validate_orientation_matrix <- function(ornt, allow_dropped = TRUE, arg_name = "ornt") {
  out <- as.matrix(ornt)
  if (ncol(out) != 2) {
    stop(sprintf("`%s` must be an n x 2 orientation matrix", arg_name))
  }

  storage.mode(out) <- "double"
  colnames(out) <- c("axis", "flip")

  axis <- out[, 1]
  flip <- out[, 2]
  axis_present <- !is.na(axis)

  if (any(axis_present & (axis %% 1 != 0))) {
    stop(sprintf("`%s[, 1]` must contain integer axis indices", arg_name))
  }

  if (any(axis_present & axis < 1)) {
    stop(sprintf("`%s[, 1]` must contain axis indices >= 1", arg_name))
  }

  if (!allow_dropped && any(!axis_present)) {
    stop(sprintf("`%s` cannot contain dropped axes (NA)", arg_name))
  }

  if (any(axis_present & is.na(flip))) {
    stop(sprintf("`%s[, 2]` cannot be NA where axis is present", arg_name))
  }

  if (any(!axis_present & !is.na(flip))) {
    stop(sprintf("`%s` uses NA axis rows; flip must also be NA for these rows", arg_name))
  }

  valid_flip <- is.na(flip) | flip %in% c(-1, 1)
  if (any(!valid_flip)) {
    stop(sprintf("`%s[, 2]` must contain only -1, 1, or NA", arg_name))
  }

  out
}

#' Normalize axis-code labels
#'
#' @param labels Label list or `NULL`.
#' @param n_axes Number of axes required.
#' @return Normalized label list.
#' @keywords internal
#' @noRd
.normalize_axcode_labels <- function(labels, n_axes) {
  if (is.null(labels)) {
    labels <- .default_axcode_labels()
  }

  if (!is.list(labels) || length(labels) < n_axes) {
    stop("`labels` must be a list with at least one (negative, positive) pair per axis")
  }

  labels <- labels[seq_len(n_axes)]
  is_valid_pair <- vapply(labels, function(x) is.character(x) && length(x) == 2, logical(1))
  if (!all(is_valid_pair)) {
    stop("Each element of `labels` must be a character vector of length 2")
  }

  label_values <- unlist(labels, use.names = FALSE)
  if (anyDuplicated(label_values)) {
    stop("Axis labels must be unique across all axes")
  }

  labels
}

#' Orientation utility functions
#'
#' Helper functions for converting between affine matrices, axis-code
#' representations, and array orientation transforms.
#'
#' Orientation matrices (`ornt`) use two columns:
#' \itemize{
#'   \item column 1 (`axis`) stores the output axis index (1-based),
#'   \item column 2 (`flip`) stores direction (`1` or `-1`).
#' }
#'
#' This is an R counterpart to NiBabel's orientation utilities, adapted to
#' `neuroim2` conventions.
#'
#' @name orientation_utils
#' @examples
#' aff <- diag(4)
#' ornt <- affine_to_orientation(aff)
#' orientation_to_axcodes(ornt)
#'
#' arr <- array(1:24, dim = c(2, 3, 4))
#' tx <- orientation_transform(
#'   axcodes_to_orientation(c("R", "A", "S")),
#'   axcodes_to_orientation(c("A", "R", "S"))
#' )
#' out <- apply_orientation(arr, tx)
#'
#' inv_aff <- orientation_inverse_affine(tx, dim(arr))
#' inv_aff
#'
#' @param affine Affine matrix of shape `(q + 1) x (p + 1)`.
#' @param tol Optional singular-value tolerance.
#' @return A `p x 2` orientation matrix with columns `axis` and `flip`.
#' @rdname orientation_utils
#' @export
affine_to_orientation <- function(affine, tol = NULL) {
  affine <- as.matrix(affine)
  if (!is.numeric(affine) || nrow(affine) < 2 || ncol(affine) < 2) {
    stop("`affine` must be a numeric matrix with at least 2 rows and 2 columns")
  }

  q <- nrow(affine) - 1
  p <- ncol(affine) - 1
  rzs <- affine[seq_len(q), seq_len(p), drop = FALSE]

  zooms <- sqrt(colSums(rzs * rzs))
  zooms[zooms == 0] <- 1
  rs <- sweep(rzs, 2, zooms, "/")

  k <- min(dim(rs))
  decomp <- svd(rs, nu = k, nv = k)

  if (is.null(tol)) {
    tol <- max(decomp$d) * max(dim(rs)) * .Machine$double.eps
  }

  keep <- decomp$d > tol
  if (any(keep)) {
    rot <- decomp$u[, keep, drop = FALSE] %*% t(decomp$v[, keep, drop = FALSE])
  } else {
    rot <- matrix(0, nrow = q, ncol = p)
  }

  out <- matrix(NA_real_, nrow = p, ncol = 2, dimnames = list(NULL, c("axis", "flip")))

  # Assign strongest-matching input axes first; ties keep original input order.
  strength <- apply(rot^2, 2, max)
  in_axes <- order(-strength, seq_along(strength))

  for (in_ax in in_axes) {
    col <- rot[, in_ax]
    if (all(abs(col) <= sqrt(.Machine$double.eps))) {
      next
    }

    out_ax <- which.max(abs(col))
    out[in_ax, 1] <- out_ax
    out[in_ax, 2] <- if (col[out_ax] < 0) -1 else 1

    # Remove the consumed output axis from further assignment.
    rot[out_ax, ] <- 0
  }

  out
}

#' Compose orientation transforms
#'
#' Returns the orientation transform that maps from `start_ornt` to `end_ornt`.
#'
#' @param start_ornt Starting orientation matrix.
#' @param end_ornt Target orientation matrix.
#' @return Orientation matrix representing `start -> end`.
#' @rdname orientation_utils
#' @export
orientation_transform <- function(start_ornt, end_ornt) {
  start_ornt <- .validate_orientation_matrix(start_ornt, allow_dropped = FALSE, arg_name = "start_ornt")
  end_ornt <- .validate_orientation_matrix(end_ornt, allow_dropped = FALSE, arg_name = "end_ornt")

  if (!identical(dim(start_ornt), dim(end_ornt))) {
    stop("`start_ornt` and `end_ornt` must have the same shape")
  }

  if (anyDuplicated(start_ornt[, 1]) || anyDuplicated(end_ornt[, 1])) {
    stop("Orientation axes must be unique")
  }

  result <- matrix(NA_real_, nrow = nrow(start_ornt), ncol = 2, dimnames = list(NULL, c("axis", "flip")))

  for (end_in_idx in seq_len(nrow(end_ornt))) {
    end_out_idx <- end_ornt[end_in_idx, 1]
    matched <- which(start_ornt[, 1] == end_out_idx)
    if (length(matched) != 1) {
      stop(sprintf("Unable to find unique axis %.0f from `end_ornt` in `start_ornt`", end_out_idx))
    }

    flip <- if (start_ornt[matched, 2] == end_ornt[end_in_idx, 2]) 1 else -1
    result[matched, ] <- c(end_in_idx, flip)
  }

  result
}

#' Flip an array over one axis
#'
#' @param arr Array.
#' @param axis Axis index (1-based).
#' @return Flipped array.
#' @keywords internal
#' @noRd
.flip_array_axis <- function(arr, axis) {
  idx <- lapply(dim(arr), seq_len)
  idx[[axis]] <- rev(idx[[axis]])
  do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
}

#' Apply an orientation matrix to array axes
#'
#' Applies flips and permutations implied by `ornt` to the first `n` dimensions
#' of `arr`, where `n = nrow(ornt)`.
#'
#' @param arr Array-like object.
#' @param ornt Orientation matrix.
#' @return Reoriented array.
#' @rdname orientation_utils
#' @export
apply_orientation <- function(arr, ornt) {
  out <- as.array(arr)
  ornt <- .validate_orientation_matrix(ornt, allow_dropped = FALSE, arg_name = "ornt")
  n <- nrow(ornt)

  if (length(dim(out)) < n) {
    stop("Data array has fewer dimensions than orientation transform")
  }

  if (anyDuplicated(ornt[, 1])) {
    stop("`ornt[, 1]` must be unique")
  }

  if (!identical(sort(as.integer(ornt[, 1])), seq_len(n))) {
    stop("`ornt[, 1]` must be a permutation of 1:n for array reorientation")
  }

  for (ax in seq_len(n)) {
    if (ornt[ax, 2] == -1) {
      out <- .flip_array_axis(out, ax)
    }
  }

  perm <- seq_along(dim(out))
  perm[seq_len(n)] <- order(ornt[, 1])
  aperm(out, perm)
}

#' Build affine that inverts an orientation transform
#'
#' Returns the affine mapping transformed array coordinates back to source array
#' coordinates for a transform encoded by `ornt`.
#'
#' @param ornt Orientation matrix.
#' @param shape Shape of source array.
#' @return Homogeneous affine matrix of size `(p + 1) x (p + 1)`.
#' @rdname orientation_utils
#' @export
orientation_inverse_affine <- function(ornt, shape) {
  ornt <- .validate_orientation_matrix(ornt, allow_dropped = FALSE, arg_name = "ornt")
  p <- nrow(ornt)

  if (anyDuplicated(ornt[, 1])) {
    stop("`ornt[, 1]` must be unique")
  }

  if (!identical(sort(as.integer(ornt[, 1])), seq_len(p))) {
    stop("`ornt[, 1]` must be a permutation of 1:p for invertible transforms")
  }

  shape <- as.numeric(shape)
  if (length(shape) < p) {
    stop("`shape` must have at least p elements, where p = nrow(ornt)")
  }

  shape <- shape[seq_len(p)]
  if (any(!is.finite(shape)) || any(shape <= 0)) {
    stop("`shape` must contain finite positive values")
  }

  axis_transpose <- as.integer(ornt[, 1])
  undo_reorder <- diag(p + 1)[c(axis_transpose, p + 1), , drop = FALSE]
  undo_flip <- diag(c(ornt[, 2], 1))

  center_trans <- -(shape - 1) / 2
  undo_flip[seq_len(p), p + 1] <- (ornt[, 2] * center_trans) - center_trans

  undo_flip %*% undo_reorder
}

#' Convert orientation matrix to axis codes
#'
#' @param ornt Orientation matrix.
#' @param labels Optional label pairs per axis.
#' @return Character vector of axis codes (positive ends), with `NA` for dropped axes.
#' @rdname orientation_utils
#' @export
orientation_to_axcodes <- function(ornt, labels = NULL) {
  ornt <- .validate_orientation_matrix(ornt, allow_dropped = TRUE, arg_name = "ornt")
  labels <- .normalize_axcode_labels(labels, nrow(ornt))

  out <- character(nrow(ornt))
  for (i in seq_len(nrow(ornt))) {
    axis <- ornt[i, 1]
    flip <- ornt[i, 2]

    if (is.na(axis)) {
      out[i] <- NA_character_
      next
    }

    axis <- as.integer(axis)
    if (axis > length(labels)) {
      stop(sprintf("Axis %d exceeds available labels", axis))
    }

    out[i] <- if (flip == 1) labels[[axis]][2] else labels[[axis]][1]
  }

  out
}

#' Convert axis codes to orientation matrix
#'
#' @param axcodes Character vector of axis codes. `NA` indicates dropped axis.
#' @param labels Optional label pairs per axis.
#' @return Orientation matrix.
#' @rdname orientation_utils
#' @export
axcodes_to_orientation <- function(axcodes, labels = NULL) {
  axcodes <- as.character(axcodes)
  labels <- .normalize_axcode_labels(labels, length(axcodes))

  out <- matrix(NA_real_, nrow = length(axcodes), ncol = 2, dimnames = list(NULL, c("axis", "flip")))

  for (i in seq_along(axcodes)) {
    code <- axcodes[i]
    if (is.na(code)) {
      next
    }

    matched <- which(vapply(labels, function(x) code %in% x, logical(1)))
    if (length(matched) != 1) {
      stop(sprintf("Axis code '%s' is invalid or ambiguous for provided labels", code))
    }

    lab <- labels[[matched]]
    out[i, ] <- c(matched, if (code == lab[1]) -1 else 1)
  }

  present_axes <- out[, 1][!is.na(out[, 1])]
  if (anyDuplicated(present_axes)) {
    stop("Axis codes map to duplicated axes, which is not a valid orientation")
  }

  out
}

#' Convert affine directly to axis codes
#'
#' @param affine Affine matrix.
#' @param labels Optional label pairs per axis.
#' @param tol Optional singular-value tolerance.
#' @return Character vector of axis codes.
#' @rdname orientation_utils
#' @export
affine_to_axcodes <- function(affine, labels = NULL, tol = NULL) {
  orientation_to_axcodes(affine_to_orientation(affine, tol = tol), labels = labels)
}

#' Get anatomical axis from signed axis index
#'
#' @param x Signed axis index in `{-3,-2,-1,1,2,3}`.
#' @return A `NamedAxis` object.
#' @keywords internal
#' @noRd
.getAxis <- function(x) {
  switch(as.character(as.integer(x)),
         "1" = LEFT_RIGHT,
         "-1" = RIGHT_LEFT,
         "2" = POST_ANT,
         "-2" = ANT_POST,
         "3" = INF_SUP,
         "-3" = SUP_INF)
}

#' Infer nearest anatomical orientation from a 3x3 transform
#'
#' @param mat33 Numeric 3x3 matrix.
#' @param tol Optional tolerance passed to `affine_to_orientation()`.
#' @return `AxisSet3D` orientation.
#' @keywords internal
#' @noRd
.nearest_anatomy_from_matrix <- function(mat33, tol = NULL) {
  mat33 <- as.matrix(mat33)
  if (!is.numeric(mat33) || !all(dim(mat33) == c(3, 3))) {
    stop("Expected a numeric 3x3 matrix")
  }

  aff <- rbind(cbind(mat33, 0), c(0, 0, 0, 1))
  ornt <- affine_to_orientation(aff, tol = tol)

  if (any(is.na(ornt[, 1]))) {
    stop("Cannot determine full anatomical orientation from rank-deficient matrix")
  }

  # findAnatomy3D() expects LPI-style axis initials.
  lpi_labels <- list(c("R", "L"), c("A", "P"), c("S", "I"))
  codes <- orientation_to_axcodes(ornt, labels = lpi_labels)
  findAnatomy3D(codes[1], codes[2], codes[3])
}

#' Find anatomical orientation from permutation matrix
#'
#' @param pmat A 3x3 permutation-like matrix.
#' @param tol Tolerance for decomposition.
#' @return An `AxisSet3D` object representing nearest anatomical orientation.
#' @keywords internal
#' @noRd
findAnatomy <- function(pmat, tol = 1e-10) {
  .nearest_anatomy_from_matrix(pmat[1:3, 1:3, drop = FALSE], tol = tol)
}

#' Find nearest valid anatomical orientation for a 4x4 transform
#'
#' @param mat44 A 4x4 transformation matrix.
#' @return An `AxisSet3D` object.
#' @keywords internal
#' @noRd
.nearestAnatomy <- function(mat44) {
  .nearest_anatomy_from_matrix(mat44[1:3, 1:3, drop = FALSE])
}
