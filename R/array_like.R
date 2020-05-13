

#' Array-like access for 4-dimensional data structures
#' @name extractor4d
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
NULL

#' Array-like access for 3-dimensional data structures
#' @name extractor3d
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param ... additional args
#' @param drop dimension
NULL



#' @export
#' @rdname extractor4d
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            ind <- grid_to_index(space(x),i)
            linear_access(x, ind)
          }
)

#' @export
#' @rdname extractor4d
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "numeric", drop="ANY"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k)) {
              k = 1:(dim(x)[3])
            }
            if (missing(m)) {
              m = 1:(dim(x)[4])
            }


            ind <- exgridToIndex4DCpp(dim(x), i,j,k,m)

            #grid <- as.matrix(expand.grid(i=i,j=j,k=k,m=m))

            ## TODO grid_to_index for 4D image doesn't work
            #ind <- .gridToIndex(dim(x), grid)

            vals <- linear_access(x,ind)
            ret <- array(vals, c(length(i), length(j), length(k), length(m)))
            if (drop) drop(ret) else ret
          }
)




#' @export
#' @rdname extractor4d
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


#' @export
#' @rdname extractor4d
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


#' @export
#' @rdname extractor4d
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


#' @export
#' @rdname extractor3d
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

#' @export
#' @rdname extractor3d
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {
            ind <- grid_to_index(x,i)
            linear_access(x, ind)
          }
)



#' @export
#' @rdname extractor3d
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

#' @export
#' @rdname extractor3d
setMethod(f="[", signature=signature(x = "ArrayLike3D", i = "missing", j = "numeric", drop="ANY"),
          def=function (x, i, j, k,  ..., drop=TRUE) {
            if (missing(k)) {
              k <- seq(1, dim(x)[3])
            }
            callGeneric(x, i=seq(1, dim(x)[1]), j, k, drop=drop, ...)
          }
)







