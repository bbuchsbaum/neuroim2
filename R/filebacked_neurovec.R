
#' @export
FileBackedNeuroVec <- function(file_name) {
  meta <- read_header(file_name)
  assert_that(length(meta@dims) == 4)
  sp <- NeuroSpace(meta@dims,meta@spacing, meta@origin,
                       meta@spatial_axes, trans(meta))

  new("FileBackedNeuroVec",
      space=sp,
      meta=meta)

}


#' @export
#' @rdname series-methods
setMethod(f="series", signature=signature(x="FileBackedNeuroVec", i="numeric"),
          def=function(x,i) {
            read_mapped_series(x@meta,i)
          })



#' @export
#' @rdname series-methods
setMethod(f="series", signature=signature(x="FileBackedNeuroVec", i="matrix"),
          def=function(x,i) {
            idx <- grid_to_index(x@space, i)
            callGeneric(x,idx)
          })


#' @export
#' @rdname sub_vector-methods
setMethod(f="sub_vector", signature=signature(x="FileBackedNeuroVec", i="numeric"),
          def=function(x, i) {
            mat <- read_mapped_vols(x@meta, i)
            sp <- add_dim(drop_dim(space(x)), length(i))
            DenseNeuroVec(mat, sp)
          })


#' as.list
#'
#' convert FileBackedNeuroVec to list of \code{\linkS4class{DenseNeuroVol}}
#' @rdname as.list-methods
#' @export
setMethod(f="as.list", signature=signature(x = "FileBackedNeuroVec"), def=function(x) {
  D4 <- dim(x)[4]

  f <- function(i) {
    drop(sub_vector(x, i))
  }

  deferred_list(lapply(seq(1, D4),
         function(i) { f } ))

})


#' extractor
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param k third index
#' @param m third index
#' @param ... additional args
#' @param drop drop dimension
setMethod(f="[", signature=signature(x = "FileBackedNeuroVec", i = "numeric", j = "missing"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
            if (missing(k) && missing(m) && nargs() == 4) {
              vals <- read_mapped_data(x@meta, i)
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

setMethod(f="[", signature=signature(x = "FileBackedNeuroVec", i = "numeric", j = "numeric"),
          def=function (x, i, j, k, m, ..., drop=TRUE) {
              if (missing(k))
                k = seq(1, dim(x)[3])
              if (missing(m)) {
                m <- seq(1, dim(x)[4])
              }

              vmat <- expand.grid(i=i, j=j, k=k, m=m)
              idx <- .gridToIndex(dim(x), vmat)
              vals <- read_mapped_data(x@meta, idx)
              ret <- array(vals, c(length(i), length(j), length(k), length(m)))
              if (drop) {
                drop(ret)
              } else {
                ret
              }

            }

)


