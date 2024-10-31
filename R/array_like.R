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

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            ind <- grid_to_index(space(x),i)
            linear_access(x, ind)
          }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "numeric", drop="ANY"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k)) {
              k = 1:(dim(x)[3])
            }
            if (missing(m)) {
              m = 1:(dim(x)[4])
            }

            ind <- exgridToIndex4DCpp(dim(x), i,j,k,m)

            vals <- linear_access(x,ind)
            ret <- array(vals, c(length(i), length(j), length(k), length(m)))
            if (drop) drop(ret) else ret
          }
)

#' @rdname extractor4d
#' @method [ ArrayLike4D
#' @export
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k) && missing(m) && nargs() == 4) {
              linear_access(x, i)
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
            callGeneric(x, i:(dim(x)[1]), j,k,m,drop=drop)
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






