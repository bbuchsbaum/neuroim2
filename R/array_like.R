#' Array-like access for 4-dimensional data structures
#'
#' @description
#' This generic function provides array-like access for 4-dimensional data structures.
#' It allows for flexible indexing and subsetting of 4D arrays or array-like objects.
#'
#' @param x The 4-dimensional object to be accessed.
#' @param i First index or dimension.
#' @param j Second index or dimension.
#' @param k Third index or dimension.
#' @param m Fourth index or dimension.
#' @param ... Additional arguments passed to methods.
#' @param drop Logical. If TRUE, the result is coerced to the lowest possible dimension.
#'
#' @return A subset of the input object, with dimensions depending on the indexing and the `drop` parameter.
#'
#' @name extractor4d
#' @rdname extractor4d
NULL

#' Array-like access for 3-dimensional data structures
#'
#' @description
#' This generic function provides array-like access for 3-dimensional data structures.
#' It allows for flexible indexing and subsetting of 3D arrays or array-like objects.
#'
#' @param x The 3-dimensional object to be accessed.
#' @param i First index or dimension.
#' @param j Second index or dimension.
#' @param k Third index or dimension.
#' @param ... Additional arguments passed to methods.
#' @param drop Logical. If TRUE, the result is coerced to the lowest possible dimension.
#'
#' @return A subset of the input object, with dimensions depending on the indexing and the `drop` parameter.
#'
#' @name extractor3d
#' @rdname extractor3d
NULL

# ------------------------------------------------------------------
# Assumed Existing Functions
# ------------------------------------------------------------------
# These functions are assumed to be correctly implemented.
# - linear_access(x, ind): Retrieves data from 'x' using linear indices 'ind'.
# - grid_to_index(space, i): Converts multi-dimensional indices to linear indices for 3D.
# - exgridToIndex4DCpp(dimensions, i, j, k, m): Converts multi-dimensional indices to linear indices for 4D.
# - space(x): Retrieves spatial information or dimension details from 'x'.

# ------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------

#' Validate Indices
#'
#' A helper function to validate numeric indices for each dimension.
#'
#' @param dimensions A numeric vector representing the dimensions of an array.
#' @param indices A list of numeric index vectors, one for each dimension to be validated.
#' @param dim_names A character vector of names corresponding to each dimension.
#'
#' @return Invisibly returns if all indices are valid. 
#'
#' @keywords internal
#' @noRd
validate_indices <- function(dimensions, indices, dim_names) {
  for (d in seq_along(indices)) {
    idx <- indices[[d]]
    if (!is.numeric(idx)) {
      stop(sprintf("Indices for dimension '%s' must be numeric.", dim_names[d]))
    }
    if (any(idx < 1 | idx > dimensions[d], na.rm = TRUE)) {
      stop(sprintf("Index out of bounds for dimension '%s'. Valid range: 1-%d, provided: %s",
                   dim_names[d], dimensions[d], paste(idx, collapse = ",")))
    }
    if (any(is.na(idx))) {
      stop(sprintf("NA values are not allowed in indices for dimension '%s'.", dim_names[d]))
    }
  }
}

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(
  f = "[",
  signature = signature(x = "ArrayLike4D", i = "matrix", j = "missing", drop = "ANY"),
  definition = function(x, i, j, k, m, ..., drop = TRUE) {
    # Error Handling
    if (!is.matrix(i)) {
      stop("When 'i' is provided as an index, it must be a matrix for 4D indexing.")
    }
    if (ncol(i) != 4) {
      stop("Matrix 'i' must have exactly 4 columns for 4D indexing (i, j, k, m).")
    }

    dims <- dim(x)
    dim_names <- c("i", "j", "k", "m")

    # Validate all columns of the matrix against corresponding dimensions
    idx_list <- list(i = i[, 1], j = i[, 2], k = i[, 3], m = i[, 4])
    validate_indices(dims, idx_list, dim_names)

    # Translate multi-dimensional indices to linear indices
    ind <- grid_to_index(space(x), i)

    # Access the data
    linear_access(x, ind)
  }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(
  f = "[",
  signature = signature(x = "ArrayLike4D", i = "numeric", j = "numeric", drop = "ANY"),
  definition = function(x, i, j, k, m, ..., drop = TRUE) {
    dims <- dim(x)
    dim_names <- c("i", "j", "k", "m")

    # Default handling for missing k and m
    if (missing(k)) {
      k <- seq_len(dims[3])
    }
    if (missing(m)) {
      m <- seq_len(dims[4])
    }

    # Validate indices
    validate_indices(dims, list(i = i, j = j, k = k, m = m), dim_names)

    # Translate multi-dimensional indices to linear indices
    ind <- exgridToIndex4DCpp(dims, i, j, k, m)

    # Access the data
    vals <- linear_access(x, ind)
    ret <- array(vals, dim = c(length(i), length(j), length(k), length(m)))

    # Apply drop logic
    if (drop) {
      return(drop(ret))
    } else {
      return(ret)
    }
  }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k) && missing(m) && nargs() == 4) {
              linear_access(x,i)
            } else {
              j <- seq(1, dim(x)[2])
              if (missing(k))
                k = seq(1, dim(x)[3])
              if (missing(m)) {
                m <- seq(1, dim(x)[4])
              }
              callGeneric(x,i,j,k,m,drop=drop)
            }
          }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "integer", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k) && missing(m) && nargs() == 4) {
              linear_access(x,i)
            } else {
              j <- seq(1, dim(x)[2])
              if (missing(k))
                k = seq(1, dim(x)[3])
              if (missing(m)) {
                m <- seq(1, dim(x)[4])
              }
              callGeneric(x,i,j,k,m,drop=drop)
            }
          }
)


#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "missing", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k)) {
              k = 1:(dim(x)[3])
            }

            if (missing(m)) {
              m = 1:(dim(x)[4])
            }

            callGeneric(x, 1:(dim(x)[1]), 1:(dim(x)[2]), k,m,drop=drop)
          }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "missing", j = "numeric"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k)) {
              k = 1:(dim(x)[3])
            }

            if (missing(m)) {
              m = 1:(dim(x)[4])
            }
            callGeneric(x, 1:(dim(x)[1]), j,k,m,drop=drop)
          }
)

#' @rdname extractor3d
#' @method [ ArrayLike3D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "numeric", j = "missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            if (missing(k) && nargs() == 4) {
              linear_access(x,i)
            } else {
              if (missing(k)) {
                k <- 1:(dim(x)[3])
              }
              callGeneric(x, i=i,  j=seq(1,dim(x)[2]), k, drop)
            }
          }
)

#' @rdname extractor3d
#' @method [ ArrayLike3D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            ind <- grid_to_index(x,i)
            linear_access(x, ind)
          }
)

#' @rdname extractor3d
#' @method [ ArrayLike3D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "missing", j = "missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            if (missing(k)) {
              idx <- seq(1, prod(dim(x)))
              callGeneric(x, idx)
            } else {
              if (missing(k)) {
                k <- seq(1, dim(x)[3])
              }
              callGeneric(x, i=seq(1, dim(x)[1]), j=seq(1, dim(x)[2]), k=k, drop=drop)
            }
          }
)

#' @rdname extractor3d
#' @method [ ArrayLike3D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "missing", j = "numeric", drop="ANY"),
          def=function (x, i, j, k,  ..., drop=TRUE) {
            if (missing(k)) {
              k <- seq(1, dim(x)[3])
            }
            callGeneric(x, i=seq(1, dim(x)[1]), j, k, drop=drop, ...)
          }
)

#' Array extraction for ClusteredNeuroVec
#'
#' @description
#' Provides array-like access to ClusteredNeuroVec objects, supporting 
#' extraction patterns like x[,,,t] to get 3D volumes at specific time points.
#'
#' @rdname extractor4d
#' @export
setMethod("[",
  signature(x = "ClusteredNeuroVec", i = "missing", j = "missing"),
  function(x, i, j, k, m, ..., drop = TRUE) {
    # Handle case where only time index is provided (x[,,,t])
    if (!missing(m) && is.numeric(m)) {
      sp3 <- dim(space(x@cvol))
      nsp <- prod(sp3)
      
      m <- as.integer(m)
      stopifnot(all(m >= 1L & m <= nrow(x@ts)))
      
      # For each selected timepoint, fill a 3D array with cluster values
      out <- lapply(m, function(ti) {
        buf <- rep.int(NA_real_, nsp)
        active <- which(x@cl_map > 0L)
        cid <- x@cl_map[active]
        buf[active] <- x@ts[ti, cid]
        array(buf, dim = sp3)
      })
      
      if (length(out) == 1 && drop) {
        out[[1]]
      } else {
        arr <- array(unlist(out, use.names = FALSE), dim = c(sp3, length(m)))
        if (drop) drop(arr) else arr
      }
    } else {
      # Delegate to generic method for other indexing patterns
      callNextMethod()
    }
  }
)
