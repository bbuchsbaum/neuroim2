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
#' @examples
#' \dontrun{
#' matchAxis("L")  # Returns LEFT_RIGHT axis
#' matchAxis("ANTERIOR")  # Returns ANT_POST axis
#' }
matchAxis <- function(firstAxis) {
  switch(toupper(firstAxis),
         "LEFT"=LEFT_RIGHT,
         "L"=LEFT_RIGHT,
         "RIGHT"=RIGHT_LEFT,
         "R"=RIGHT_LEFT,
         "ANTERIOR"=ANT_POST,
         "A"=ANT_POST,
         "POSTERIOR"=POST_ANT,
         "P"=POST_ANT,
         "INFERIOR"=INF_SUP,
         "I"=INF_SUP,
         "SUPERIOR"=SUP_INF,
         "S"=SUP_INF)

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
AxisSet1D <- function(i) {
  new("AxisSet1D", ndim=as.integer(1), i=i)
}

#' Create a two-dimensional axis set
#'
#' @param i A NamedAxis object representing the first axis
#' @param j A NamedAxis object representing the second axis
#' @return An AxisSet2D object
#' @keywords internal
AxisSet2D <- function(i, j) {
	new("AxisSet2D", ndim=as.integer(2), i=i, j=j)
}

#' Create a three-dimensional axis set
#'
#' @param i A NamedAxis object representing the first axis
#' @param j A NamedAxis object representing the second axis
#' @param k A NamedAxis object representing the third axis
#' @return An AxisSet3D object
#' @keywords internal
AxisSet3D <- function(i, j, k) {
	new("AxisSet3D", ndim=as.integer(3), i=i, j=j, k=k)
}

#' Get permutation matrix from axis set
#'
#' @param x An AxisSet2D object
#' @param ... Additional arguments (not used)
#' @return A matrix representing the axis directions
#' @export
#' @examples
#' \dontrun{
#' axes <- AxisSet2D(LEFT_RIGHT, ANT_POST)
#' perm_mat(axes)
#' }
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
#' @examples
#' \dontrun{
#' axes <- AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP)
#' perm_mat(axes)
#' }
setMethod(f="perm_mat", signature=signature(x = "AxisSet3D"),
          def=function(x, ...) {
            cbind(x@i@direction, x@j@direction, x@k@direction)
          })

#' Get permutation matrix from axis set
#'
#' @param x A NeuroSpace object
#' @param ... Additional arguments (not used)
#' @return A matrix representing the axis directions
#' @export
#' @examples
#' \dontrun{
#' ns <- NeuroSpace(axes=AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP))
#' perm_mat(ns)
#' }
setMethod(f="perm_mat", signature=signature(x = "NeuroSpace"),
          def=function(x, ...) {
            callGeneric(x@axes)
          })

#' Drop dimension from axis set
#'
#' @param x An AxisSet2D or AxisSet3D object
#' @param dimnum Numeric index of dimension to drop
#' @return An axis set with one fewer dimension
#' @export
#' @examples
#' \dontrun{
#' axes3d <- AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP)
#' axes2d <- drop_dim(axes3d, 1)  # Drop first dimension
#' }
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

#' Drop dimension from axis set
#'
#' @param x An AxisSet2D object
#' @param dimnum Numeric index of dimension to drop (optional)
#' @return An axis set with one fewer dimension
#' @export
#' @examples
#' \dontrun{
#' axes2d <- AxisSet2D(LEFT_RIGHT, ANT_POST)
#' axes1d <- drop_dim(axes2d)  # Drop first dimension
#' }
setMethod(f="drop_dim", signature=signature(x = "AxisSet2D", dimnum="missing"),
          def=function(x, dimnum) {
            AxisSet1D(x@i)
          })

#' Drop dimension from axis set
#'
#' @param x An AxisSet3D object
#' @param dimnum Numeric index of dimension to drop
#' @return An axis set with one fewer dimension
#' @export
#' @examples
#' \dontrun{
#' axes3d <- AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP)
#' axes2d <- drop_dim(axes3d, 1)  # Drop first dimension
#' }
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
#' @return An axis set with one fewer dimension
#' @export
#' @examples
#' \dontrun{
#' axes3d <- AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP)
#' axes2d <- drop_dim(axes3d)  # Drop first dimension
#' }
setMethod(f="drop_dim", signature=signature(x = "AxisSet3D", dimnum="missing"),
           def=function(x, dimnum) {
             AxisSet2D(x@i, x@j)
           })

#' Get number of dimensions in axis set
#'
#' @param x An AxisSet object
#' @param ... Additional arguments (not used)
#' @return Integer number of dimensions
#' @export
#' @examples
#' \dontrun{
#' axes <- AxisSet2D(LEFT_RIGHT, ANT_POST)
#' ndim(axes)  # Returns 2
#' }
setMethod(f="ndim",signature(x= "AxisSet"), function(x, ...) { x@ndim })

#' Print method for NamedAxis objects
#'
#' @param x A NamedAxis object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis
#' @keywords internal
setMethod(f="print_", signature=signature("NamedAxis"),
    def=function(x, ...) {
        x@axis
    })

#' Show method for NamedAxis objects
#'
#' @param object A NamedAxis object
#' @export
setMethod(f="show", signature("NamedAxis"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("NamedAxis"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("─", 30), collapse="")), "\n")
        cat(crayon::white(object@axis), "\n")
    })

#' Show method for AxisSet1D objects
#'
#' @param object An AxisSet1D object
#' @export
setMethod(f="show", signature=signature("AxisSet1D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet1D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("─", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis:"), crayon::white(object@i@axis), "\n")
    })

#' Print method for AxisSet2D objects
#'
#' @param x An AxisSet2D object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis set
#' @keywords internal
setMethod(f="print_", signature=signature("AxisSet2D"),
    def=function(x, ...) {
        paste(x@i@axis, "×", x@j@axis)
    })

#' Show method for AxisSet2D objects
#'
#' @param object An AxisSet2D object
#' @export
setMethod(f="show", signature=signature("AxisSet2D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet2D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("─", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis 1:"), crayon::white(object@i@axis), "\n")
        cat(crayon::yellow("Axis 2:"), crayon::white(object@j@axis), "\n")
    })

#' Print method for AxisSet3D objects
#'
#' @param x An AxisSet3D object
#' @param ... Additional arguments (not used)
#' @return Character string representing the axis set
#' @keywords internal
setMethod(f="print_", signature=signature("AxisSet3D"),
    def=function(x, ...) {
        paste(x@i@axis, "×", x@j@axis, "×", x@k@axis)
    })

#' Show method for AxisSet3D objects
#'
#' @param object An AxisSet3D object
#' @export
setMethod(f="show", signature=signature("AxisSet3D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet3D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("─", 30), collapse="")), "\n")
        cat(crayon::yellow("Axis 1:"), crayon::white(object@i@axis), "\n")
        cat(crayon::yellow("Axis 2:"), crayon::white(object@j@axis), "\n")
        cat(crayon::yellow("Axis 3:"), crayon::white(object@k@axis), "\n")
    })

#' Show method for AxisSet4D objects
#'
#' @param object An AxisSet4D object
#' @export
setMethod(f="show", signature=signature("AxisSet4D"),
    def=function(object) {
        header <- crayon::bold(crayon::blue("AxisSet4D"))
        cat(header, "\n")
        cat(crayon::silver(paste(rep("─", 30), collapse="")), "\n")
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
#' \dontrun{
#' # Create orientation with default LPI axes
#' orient <- findAnatomy3D()
#' # Create orientation with custom axes
#' orient <- findAnatomy3D("R", "A", "S")
#' }
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
matchAnatomy3D <- function(axis1, axis2, axis3) {
	for (orient in OrientationList3D) {
		if (identical(orient@i,axis1) && identical(orient@j,axis2) && identical(orient@k, axis3)) {
			return(orient)
		}
	}

	stop("did not find legal matching anatomical orientation for axes: ",
			axis1, axis2, axis3)

}

#' Find anatomical orientation from permutation matrix
#'
#' Determines the anatomical orientation corresponding to a given permutation matrix.
#'
#' @param pmat A 3x3 permutation matrix
#' @param tol Tolerance for numerical comparisons (default: 1e-10)
#' @return An AxisSet3D object representing the anatomical orientation
#' @export
#' @examples
#' \dontrun{
#' # Create a permutation matrix
#' pmat <- matrix(c(1,0,0, 0,1,0, 0,0,1), 3, 3)
#' # Find corresponding anatomical orientation
#' orient <- findAnatomy(pmat)
#' }
findAnatomy <- function(pmat, tol=1e-10) {
  mat33 <- pmat[1:3, 1:3]
  #mat33 <- sweep(mat33, 2, sqrt(apply(mat33 * mat33, 2, sum)), "/")
  icol <- mat33[,1]
  jcol <- mat33[,2]
  kcol <- mat33[,3]

  # if (!sum(mat33 == 1) == 3) {
  #   mat33 = apply(mat33, 2, function(x) (x+runif(3)*.1))
  #   apply(mat33, 2, function(x) x/sum(x))
  # }

  ## normalize icol
  icol <- icol / sqrt(sum(icol^2))

  ## normalize jcol
  jcol <- jcol / sqrt(sum(jcol^2))

  #' Orthogonalize two vectors
  #'
  #' Internal helper function to make two vectors orthogonal while preserving
  #' the direction of the first vector.
  #'
  #' @param col1 First vector to preserve
  #' @param col2 Second vector to orthogonalize
  #' @return Orthogonalized second vector
  #' @keywords internal
  #' @noRd
  orthogonalize <- function(col1, col2) {
    dotp <- sum(col1*col2)
    #print(dotp)
    if (abs(dotp) > 1.e-4) {
      col2 <- col2 - (dotp * col1)
      norm <- sqrt(sum(col2^2))
      col2 <- col2 / norm
    }
    col2
  }

  jcol <- orthogonalize(icol, jcol)

  knorm <- sqrt(sum(kcol^2))
  if (knorm == 0.0) {
    kcol[1] <- icol[2] * jcol[3] - icol[3] * jcol[2]
    kcol[2] <- icol[3] * jcol[1] - jcol[3] * icol[1]
    kcol[3] <- icol[1] * jcol[2] - icol[2] * jcol[1]
  } else {
    kcol <- kcol / knorm
  }

  ## orthogonalize k to i
  kcol <- orthogonalize(icol, kcol)
  kcol <- orthogonalize(jcol, kcol)

  Q <- cbind(icol, jcol, kcol)
  P <- matrix(0,3,3)
  detQ <- det(Q)
  if (detQ == 0.0) {
    stop("invalid matrix input, determinant is 0")
  }


  vbest = -666
  ibest = 1
  pbest = 1
  qbest = 1
  rbest = 1

  jbest = 2
  kbest = 3
  for (i in 1:3) {
    for (j in 1:3) {
      if (i == j) next
      for (k in 1:3) {
        if (i == k || j == k) next
          P <- matrix(0,3,3)
          for (p in c(-1,1)) {
            for (q in c(-1,1)) {
              for (r in c(-1,1)) {
                P[1, i] <- p
                P[2, j] <- q
                P[3, k] <- r
                detP <- det(P)
                if (detP * detQ <= 0.0) next
                M <- P %*% Q
                crit <- sum(diag(M))
                if (crit > vbest) {
                  vbest = crit
                  ibest = i
                  jbest = j
                  kbest = k
                  pbest = p
                  qbest = q
                  rbest = r
                }
              }
            }
          }
      }
    }
  }

  .getAxis <- function(num) {
    switch(as.character(as.integer(num)),
           "1"=LEFT_RIGHT,
           "-1"=RIGHT_LEFT,
           "2"=POST_ANT,
           "-2"=ANT_POST,
           "3"=INF_SUP,
           "-3"=SUP_INF)
  }

  ax1 <- .getAxis(ibest*pbest)
  ax2 <- .getAxis(jbest*qbest)
  ax3 <- .getAxis(kbest*rbest)

  AxisSet3D(ax1,ax2,ax3)
}

#' Get anatomical axis from direction vector
#'
#' Internal helper function to convert a direction vector to an anatomical axis.
#'
#' @param x Numeric vector representing direction
#' @return A NamedAxis object
#' @keywords internal
.getAxis <- function(x) {
  switch(as.character(as.integer(x)),
         "1"=LEFT_RIGHT,
         "-1"=RIGHT_LEFT,
         "2"=POST_ANT,
         "-2"=ANT_POST,
         "3"=INF_SUP,
         "-3"=SUP_INF)
}

#' Find nearest valid anatomical orientation for a transformation matrix
#'
#' Internal helper function that finds the closest valid anatomical orientation
#' for a given 4x4 transformation matrix.
#'
#' @param mat44 A 4x4 transformation matrix
#' @return An AxisSet3D object representing the nearest valid anatomical orientation
#' @noRd
#' @keywords internal
.nearestAnatomy <- function(mat44) {
  mat33 <- mat44[1:3, 1:3]
  #mat33 <- sweep(mat33, 2, sqrt(apply(mat33 * mat33, 2, sum)), "/")
  icol <- mat33[,1]
  jcol <- mat33[,2]
  kcol <- mat33[,3]

  # if (!sum(mat33 == 1) == 3) {
  #   mat33 = apply(mat33, 2, function(x) (x+runif(3)*.1))
  #   apply(mat33, 2, function(x) x/sum(x))
  # }

  ## normalize icol
  icol <- icol / sqrt(sum(icol^2))

  ## normalize jcol
  jcol <- jcol / sqrt(sum(jcol^2))

  #' Orthogonalize two vectors
  #'
  #' Internal helper function to make two vectors orthogonal while preserving
  #' the direction of the first vector.
  #'
  #' @param col1 First vector to preserve
  #' @param col2 Second vector to orthogonalize
  #' @return Orthogonalized second vector
  #' @keywords internal
  #' @noRd
  orthogonalize <- function(col1, col2) {
    dotp <- sum(col1*col2)
    #print(dotp)
    if (abs(dotp) > 1.e-4) {
      col2 <- col2 - (dotp * col1)
      norm <- sqrt(sum(col2^2))
      col2 <- col2 / norm
    }
    col2
  }

  jcol <- orthogonalize(icol, jcol)

  knorm <- sqrt(sum(kcol^2))
  if (knorm == 0.0) {
    kcol[1] <- icol[2] * jcol[3] - icol[3] * jcol[2]
    kcol[2] <- icol[3] * jcol[1] - jcol[3] * icol[1]
    kcol[3] <- icol[1] * jcol[2] - icol[2] * jcol[1]
  } else {
    kcol <- kcol / knorm
  }

  ## orthogonalize k to i
  kcol <- orthogonalize(icol, kcol)
  kcol <- orthogonalize(jcol, kcol)

  Q <- cbind(icol, jcol, kcol)
  P <- matrix(0,3,3)
  detQ <- det(Q)
  if (detQ == 0.0) {
    stop("invalid matrix input, determinant is 0")
  }

  vbest = -666
  ibest = 1
  pbest = 1
  qbest = 1
  rbest = 1

  jbest = 2
  kbest = 3
  for (i in 1:3) {
    for (j in 1:3) {
      if (i == j) next
      for (k in 1:3) {
        if (i == k || j == k) next
          P <- matrix(0,3,3)
          for (p in c(-1,1)) {
            for (q in c(-1,1)) {
              for (r in c(-1,1)) {
                P[1, i] <- p
                P[2, j] <- q
                P[3, k] <- r
                detP <- det(P)
                if (detP * detQ <= 0.0) next
                M <- P %*% Q
                crit <- sum(diag(M))
                if (crit > vbest) {
                  vbest = crit
                  ibest = i
                  jbest = j
                  kbest = k
                  pbest = p
                  qbest = q
                  rbest = r
                }
              }
            }
          }
      }
    }
  }

  .getAxis <- function(num) {
    switch(as.character(as.integer(num)),
           "1"=LEFT_RIGHT,
           "-1"=RIGHT_LEFT,
           "2"=POST_ANT,
           "-2"=ANT_POST,
           "3"=INF_SUP,
           "-3"=SUP_INF)
  }

  ax1 <- .getAxis(ibest*pbest)
  ax2 <- .getAxis(jbest*qbest)
  ax3 <- .getAxis(kbest*rbest)

  AxisSet3D(ax1,ax2,ax3)
}
