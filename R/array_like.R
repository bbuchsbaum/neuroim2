#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m third index
#' @param ... additional args
#' @param drop drop dimension
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k) && missing(m) && nargs() == 4) {
              vals <- linear_access(x@meta, i)
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


#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "ArrayLike4D", i = "numeric", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k)) {
              k = 1:(dim(x)[3])
            }
            if (missing(m)) {
              m = 1:(dim(x)[4])
            }

            callGeneric(x, i, 1:(dim(x)[2]),k,m,drop=drop)
          }
)


#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
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


#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
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


