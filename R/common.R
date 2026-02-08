#' @export
#' @rdname read_columns-methods
setMethod(f="read_columns", signature=c(x="ColumnReader", column_indices="integer"),
          def=function(x, column_indices) {
            x@reader(column_indices)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="integer", FUN="function"),
          def=function(x, fac, FUN) {
            callGeneric(x,as.factor(fac), FUN)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="factor", FUN="missing"),
          def=function(x, fac) {
            split_reduce(x, fac, mean)
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "matrix", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            if (length(fac) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }

            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(names(ind), function(lev) {
              future.apply::future_apply(x[ind[[lev]],,drop=FALSE], 2, FUN)
            }))

            row.names(out) <- levels(fac)
            out
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "NeuroVec", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            if (length(fac) == prod(dim(x)[1:3])) {
              split_by_voxel <- TRUE
            } else if (length(fac) == dim(x)[4]) {
              split_by_row <- TRUE
            } else {
              stop(paste("length of 'fac' must be equal to number of voxels or to number of volumes"))
            }

            if (split_by_voxel) {

              ind <- split(seq_along(fac), fac)
              out <- do.call(rbind, lapply(names(ind), function(lev) {
                #idx <- which(fac == lev)
                mat <- series(x, ind[[lev]])
                apply(mat, 1, FUN)
              }))

              row.names(out) <- levels(fac)
              out
            } else {
              m <- as.matrix(x)
              callGeneric(m, fac, FUN)
            }
          })


#' @export
#' @rdname split_reduce-methods
setMethod(f="split_reduce", signature=signature(x = "NeuroVec", fac="factor", FUN="missing"),
          def=function(x, fac, FUN) {
            if (length(fac) == prod(dim(x)[1:3])) {
              split_by_voxel <- TRUE
            } else if (length(fac) == dim(x)[4]) {
              split_by_row <- TRUE
            } else {
              stop(paste("length of 'fac' must be equal to number of voxels or to number of volumes"))
            }

            if (split_by_voxel) {

              ind <- split(seq_along(fac), fac)
              out <- do.call(rbind, lapply(names(ind), function(lev) {
                rowMeans(series(x, ind[[lev]]))
              }))

              row.names(out) <- levels(fac)
              out
            } else {
              m <- as.matrix(x)
              callGeneric(m, fac)
            }
          })




#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "matrix", f="factor", center="logical", scale="logical"),
          def=function(x, f, center=TRUE, scale=TRUE) {
            if (length(f) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }

            out <- matrix(0, nrow(x), ncol(x))
            ind <- split(seq_along(f), f)

            for (lev in names(ind)) {
              keep <- ind[[lev]]
              xs <- base::scale(x[keep,,drop=FALSE], center=center, scale=scale)
              out[keep,] <- xs
            }

            out
          })



#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "matrix", f="factor", center="missing", scale="missing"),
          def=function(x, f) {
            callGeneric(x,f, TRUE, TRUE)
          })


#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="logical", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, center, TRUE)
          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="missing", scale="missing"),
          def=function(x, f) {
            callGeneric(x, f, TRUE, TRUE)

          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="logical", scale="missing"),
          def=function(x, f, center) {
            callGeneric(x, f, center, TRUE)

          })


#' @export
#' @rdname split_scale-methods
setMethod(f="split_scale", signature=signature(x = "DenseNeuroVec", f="factor", center="logical", scale="logical"),
          def=function(x, f, center, scale) {
            m <- callGeneric(t(as.matrix(x)), f, center, scale)
            NeuroVec(m, space(x))
          })



#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="DenseNeuroVec", center="logical", scale="logical"),
          def=function(x, center, scale) {
            d <- dim(x)
            nv <- prod(d[1:3])
            nt <- d[4]
            # Reshape directly to voxels x time — no transpose needed
            M <- matrix(x@.Data, nrow = nv, ncol = nt)
            if (center) {
              M <- M - rowMeans(M)
            }
            if (scale) {
              rsd <- sqrt(rowSums(M * M) / (nt - 1L))
              rsd[rsd == 0] <- 1
              M <- M / rsd
            }
            dim(M) <- d
            new("DenseNeuroVec", .Data = M, space = space(x))
          })

#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="logical", scale="logical"),
          def=function(x, center, scale) {
            M <- as.matrix(x)
            Ms <- base::scale(t(M), center, scale)
            NeuroVec(Ms, space(x))

          })


#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="missing", scale="logical"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, scale)
          })



#' @export
#' @rdname scale_series-methods
setMethod(f="scale_series", signature=signature(x="NeuroVec", center="missing", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, TRUE)
          })


#' .isExtension
#' @keywords internal
#' @noRd
.isExtension <- function(fname, extension) {
  last <- substr(fname, nchar(fname)+1 - nchar(extension), nchar(fname))
  return(last==extension)
}

#' .concat4D
#' @keywords internal
#' @noRd
.concat4D <- function(x, y, ...) {
  rest <- list(...)

  D <- dim(x)[1:3]

  rvols <- lapply(rest, function(z) {
    stopifnot(length(dim(z)) >= 3)
    stopifnot(identical(D, dim(z)[1:3]))
    as(z, "matrix")
  })

  clist <- if (length(rvols) > 0) {
    c(list(as(x, "matrix"), as(y, "matrix")), rvols)
  } else {
    list(as(x, "matrix"), as(y, "matrix"))
  }


  out <- do.call(cbind, clist)

  nvols <- ncol(out)
  new.dim <- c(D, nvols)

  nspace <-
    NeuroSpace(
      new.dim,
      origin = origin(x@space),
      spacing = spacing(x@space),
      axes = axes(x@space),
      trans = trans(x@space)
    )

  DenseNeuroVec(out, nspace)

}




#' .gridToIndex3D
#' @importFrom assertthat assert_that
#' @keywords internal
#' @noRd
.gridToIndex3D <- function(dimensions, voxmat) {
	assert_that(length(dimensions) == 3)
  if (is.vector(voxmat)) {
    assert_that(length(voxmat) == 3)
    voxmat <- matrix(voxmat, 1,3)
  }

  assert_that(ncol(voxmat) == 3)
  gridToIndex3DCpp(dimensions, voxmat)

}

#' .gridToIndex
#' @keywords internal
#' @importFrom purrr map_dbl
#' @noRd
.gridToIndex <- function(dimensions, vmat) {
  vmat <- as.matrix(vmat)
  assert_that(length(dimensions) == ncol(vmat), msg=paste("length(dimensions) not equal to ncol(vmat): ", length(dimensions), "!=", ncol(vmat)))

	# D <- Reduce("*", dimensions, accumulate=TRUE)
	# apply(vmat, 1, function(vox) {
	# 	sum(map_dbl(length(D):2, function(i) {
	# 		D[i-1]*(vox[i]-1)
	# 	})) + vox[1]
	# })

  gridToIndexCpp(as.integer(dimensions), vmat)

}

#' .indexToGrid
#' @keywords internal
#' @noRd
.indexToGrid <- function(idx, array.dim) {
  assert_that(all(idx > 0 & idx <= prod(array.dim)))
  assert_that(length(array.dim) <= 5)
  indexToGridCpp(idx, array.dim)

}



#' .getRStorage
#' @keywords internal
#' @noRd
.getRStorage <- function(data_type) {
  dtype_upper <- toupper(data_type)
  if (any(dtype_upper == c("BINARY", "BYTE", "UBYTE", "SHORT", "INTEGER", "INT", "LONG"))) {
    "integer"
  } else if (any(dtype_upper == c("FLOAT", "DOUBLE"))) {
    "double"
  } else {
	  stop(paste("unrecognized data type", data_type))
  }
}


#' @noRd
.isSigned <- function(data_type) {
  if (data_type == "UBYTE") {
    FALSE
  } else {
    TRUE
  }
}



#' @keywords internal
#' @importFrom mmap int8 uint8 int16 int32 real32 real64
#' @importFrom mmap mmap char mmapFlags munmap
#' @noRd
.getMMapMode <- function(code) {
	if (code == "UNKNOWN") {
		stop(paste(".getMMapMode: no memory map mode for UNKNOWN data type: ", code))
	} else if (code == "BINARY") {
		mmap::int8()
	} else if (code == "UBYTE") {
	  mmap::uint8()
	} else if(code == "SHORT") {
	  mmap::int16()
	} else if(code == "INT") {
	  mmap::int32()
	} else if (code == "FLOAT") {
	  mmap::real32()
	} else if (code == "DOUBLE") {
	  mmap::real64()
	} else {
		stop(paste(".getMMapMode: unsupported data type: ", code))
	}
}


#' .getDataStorage
#' @keywords internal
#' @noRd
.getDataStorage <- function(code) {
  if (code == 0) {
    return("UNKNOWN")
  } else if (code == 1) {
    return("BINARY")
  } else if (code == 2) {
    return("UBYTE")
  } else if(code == 4) {
    return("SHORT")
  } else if(code == 8) {
    return("INT")
  } else if (code == 16) {
    return("FLOAT")
  } else if (code == 64) {
    return("DOUBLE")
  } else {
    stop(paste("nifti(getDataStorage): unsupported data type: ", code))
  }
}

#' .getDataCode
#' @keywords internal
#' @noRd
.getDataCode <- function(data_type) {
  if (data_type == "UNKNOWN") {
    return(0)
  }else if (data_type == "BINARY") {
    return(1)
  } else if (data_type == "UBYTE") {
    return(2)
  } else if(data_type == "SHORT") {
    return(4)
  } else if(data_type == "INT") {
    return(8)
  } else if (data_type == "FLOAT") {
    return(16)
  } else if (data_type == "DOUBLE") {
    return(64)
  } else {
    stop(paste("getDataCode: unsupported data type: ", data_type))
  }
}

#' .getDataSize
#' @keywords internal
#' @noRd
.getDataSize <- function(data_type) {
  if (data_type == "BINARY") {
    return(1)
  } else if (data_type == "BYTE") {
	  return(1)
  } else if (data_type == "UBYTE") {
    return(1)
  } else if (data_type == "SHORT") {
    return(2)
  } else if (data_type == "INTEGER") {
    return(4)
  } else if (data_type == "INT") {
    return(4)
  } else if (data_type == "FLOAT") {
    return(4)
  } else if (data_type == "DOUBLE") {
    return(8)
  } else if (data_type == "LONG") {
    return(8)
  }

  stop(paste("unrecognized data type: ", data_type))
}

#' .getEndian
#' @keywords internal
#' @noRd
.getEndian <- function(conn) {
  #try little endian
  endian <- "little"

  hsize <- readBin(conn, integer(), 1, endian=endian)
  if (hsize != 348) {
    # might be bigendian
    endian <- "big"
    seek(conn, 0)
    hsize <- readBin(conn, integer(), 1, endian=endian)
    if (hsize != 348) {
      stop("nifti(getEndian): header size is not 348, invalid header.")
    }
  }

  return(endian)
}



#' @keywords internal
#' @noRd
.niftiExt <- function(filetype) {

  extensions <- list()

  if (filetype == "nifti-single") {
    extensions[["header"]]  <- "nii"
    extensions[["data"]] <- "nii"
  }
  else if (filetype == "nifti-pair") {
    extensions[["header"]]  <- "hdr"
    extensions[["data"]] <- "img"
  }
  else if (filetype == "nifti-gz") {
    extensions[["header"]]  <- "nii.gz"
    extensions[["data"]] <- "nii.gz"
  } else {
    stop(paste("unsupported filetype: ", filetype))
  }

  return(extensions)
}

#' Convert a Transformation Matrix to a Quaternion Representation
#'
#' @description
#' Extracts the rotation and scaling components from a 3x3 (or 4x4) transformation
#' matrix, normalizes them, and computes the corresponding quaternion parameters
#' and a sign factor (`qfac`) indicating whether the determinant is negative.
#'
#' @details
#' This function first checks and corrects for zero-length axes in the upper-left
#' corner of the matrix, then normalizes each column to extract the pure rotation.
#' If the determinant of the rotation submatrix is negative, the \code{qfac} is set
#' to \code{-1}, and the third column is negated. Finally, the quaternion parameters
#' (\eqn{a, b, c, d}) are computed following standard NIfTI-1 conventions for
#' representing the rotation in 3D.
#'
#' @param mat A numeric matrix with at least the top-left 3x3 portion containing
#'   rotation/scaling. Often a 4x4 affine transform, but only the 3x3 top-left
#'   submatrix is used in practice.
#'
#' @return A named \code{list} with two elements:
#'   \describe{
#'     \item{\code{quaternion}}{A numeric vector of length 3, \eqn{(b, c, d)},
#'       which—together with \eqn{a} derived internally—represents the rotation.}
#'     \item{\code{qfac}}{Either \code{+1} or \code{-1}, indicating whether the
#'       determinant of the rotation submatrix is positive or negative, respectively.}
#'   }
#'
#' @seealso
#' \code{\link{quaternToMatrix}} for the inverse operation, converting
#' quaternion parameters back to a transform matrix.
#'
#' @references
#' - Cox RW. *Analysis of Functional NeuroImages* (AFNI) and NIfTI-1 quaternion
#'   conventions. \url{https://afni.nimh.nih.gov}
#'
#' @export
matrixToQuatern <- function(mat) {
  xd <- sqrt(drop(crossprod(mat[1:3,1])))
  yd <- sqrt(drop(crossprod(mat[1:3,2])))
  zd <- sqrt(drop(crossprod(mat[1:3,3])))

  if (xd == 0) { mat[1,1] = 1; mat[2:3,1] = 0; xd = 1; }
  if (yd == 0) { mat[2,2] = 1; mat[c(1,3),2] = 0; yd = 1; }
  if (zd == 0) { mat[3,3] = 1; mat[1:2,3] = 0; zd = 1; }

  rmat = mat[1:3, 1:3]
  rmat[,1] = rmat[,1]/xd
  rmat[,2] = rmat[,2]/yd
  rmat[,3] = rmat[,3]/zd

  ####### skipping orthogonalization of columns

  #################################################

  zd = det(rmat)
  qfac = 1

  if (zd > 0) {
    qfac = 1
  } else {
    qfac = -1
    rmat[1:3,3] = -rmat[1:3,3]
  }

  # compute quaternion parameters

  a = rmat[1,1] + rmat[2,2] + rmat[3,3] + 1

  if (a > .5) {
    a = .5 * sqrt(a)
    b = 0.25 * (rmat[3,2]-rmat[2,3]) / a
    c = 0.25 * (rmat[1,3]-rmat[3,1]) / a
    d = 0.25 * (rmat[2,1]-rmat[1,2]) / a
   } else {
     xd = 1.0 + rmat[1,1] - (rmat[2,2]+rmat[3,3])
     yd = 1.0 + rmat[2,2] - (rmat[1,1]+rmat[3,3])
     zd = 1.0 + rmat[3,3] - (rmat[1,1]+rmat[2,2])
     if( xd > 1.0 ){
       b = 0.5 * sqrt(xd)
       c = 0.25* (rmat[1,2]+rmat[2,1]) / b
       d = 0.25* (rmat[1,3]+rmat[3,1]) / b
       a = 0.25* (rmat[3,2]-rmat[2,3]) / b
     } else if( yd > 1.0 ){
       c = 0.5 * sqrt(yd) ;
       b = 0.25* (rmat[1,2]+rmat[2,1]) / c
       d = 0.25* (rmat[2,3]+rmat[3,2]) / c
       a = 0.25* (rmat[1,3]-rmat[3,1]) / c
     } else {
       d = 0.5 * sqrt(zd) ;
       b = 0.25* (rmat[1,3]+rmat[3,1]) / d
       c = 0.25* (rmat[2,3]+rmat[3,2]) / d
       a = 0.25* (rmat[2,1]-rmat[1,2]) / d
     }
     if( a < 0.0 ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

  return(list(quaternion=c(b,c,d), qfac=qfac))

}


#' Convert Quaternion Parameters to a Transformation Matrix
#'
#' @description
#' Given a quaternion \code{(b, c, d)}, a scalar offset (origin), voxel step sizes,
#' and the \code{qfac} sign, reconstructs a 4x4 affine matrix representing rotation,
#' scaling, and translation as used in NIfTI-1.
#'
#' @details
#' This function uses the quaternion formalism common in neuroimaging, adding the
#' offset (translation) into the 4th column, and applying the voxel sizes along
#' each axis. If \code{qfac} is \code{-1}, the \eqn{z} scale is negated. The
#' resulting 4x4 matrix is typically used as an affine transform for voxel-to-world
#' coordinate mapping.
#'
#' @param quat A numeric vector of length 3 containing the quaternion parameters
#'   \eqn{(b, c, d)}. The scalar part \eqn{a} is computed internally.
#' @param origin A numeric vector of length 3 specifying the translation components
#'   (often the real-space origin or offset).
#' @param stepSize A numeric vector of length 3 giving the voxel dimensions along
#'   each axis (e.g., \code{(dx, dy, dz)}).
#' @param qfac Either \code{+1} or \code{-1}, indicating the sign from the
#'   determinant check in \code{\link{matrixToQuatern}}.
#'
#' @return A 4x4 numeric affine transformation matrix. The top-left 3x3 submatrix
#'   encodes rotation and scaling, and the 4th column encodes translation.
#'
#' @seealso
#' \code{\link{matrixToQuatern}} for converting a matrix back to quaternion form.
#'
#' @export
quaternToMatrix <- function(quat, origin, stepSize, qfac) {
  mat <- matrix(0, 4,4)
  mat[4,] <- c(0,0,0,1)

  a <- 1 - sum(quat^2)
  if (a < 1e-07) {
    a <- 1 /(sqrt(sum(quat^2)))
    quat <- quat*a
    a <- 0
  } else {
    a <- sqrt(a)
  }

  stepSize <- ifelse(stepSize > 0, stepSize, 1)
  xd <- stepSize[1]
  yd <- stepSize[2]
  zd <- stepSize[3]

  if (qfac < 0) {
    zd <- -zd
  }

  b <- quat[1]
  c <- quat[2]
  d <- quat[3]


  mat[1,1] <- (a*a+b*b-c*c-d*d) * xd
  mat[1,2] <- 2 * (b*c-a*d)     * yd
  mat[1,3] <- 2 * (b*d+a*c)     * zd
  mat[2,1] <- 2 * (b*c+a*d)     * xd
  mat[2,2] <- (a*a+c*c-b*b-d*d) * yd
  mat[2,3] <- 2 * (c*d-a*b)     * zd
  mat[3,1] <- 2 * (b*d-a*c)     * xd
  mat[3,2] <- 2 * (c*d+a*b)     * yd
  mat[3,3] <- (a*a+d*d-c*c-b*b) * zd

  mat[1:3,4] <- origin

  return(mat)
}


# @rdname internal-methods
# @keywords internal
# .makeMMap <- function(meta) {
#   nels <- prod(meta@Dim[1:4])
#
#   if (.Platform$endian != meta@endian) {
#     ## read raw bytes
#     rawbytes <- mmap::mmap(meta@data_file, mode=mmap::char(), prot=mmap::mmapFlags("PROT_READ"))
#     rawbytes <- rawbytes[(meta@data_offset+1):length(rawbytes)]
#
#     mmap::munmap(rawbytes)
#     readBin(rawbytes, what=.getRStorage(meta@data_type), size=.getDataSize(meta@data_type), n=nels, endian=meta@endian)
#   } else {
#     #mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), off=meta@data_offset,prot=mmap::mmapFlags("PROT_READ"),flags=mmap::mmapFlags("MAP_PRIVATE"))
#     ret <- mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), prot=mmap::mmapFlags("PROT_READ"))
#     offset <- meta@data_offset/.getDataSize(meta@data_type) + 1
#     vals <- ret[offset:nels]
#     mmap::munmap(ret)
#     vals
#   }
#
#
# }
