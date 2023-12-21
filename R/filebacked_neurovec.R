
#' FileBackedNeuroVec
#'
#' Construct a \code{FileBackedNeuroVec} instance
#'
#' @param file_name the name of the image file
#' @export
#' @return a new instance of type \code{FileBackedNeuroVec}
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
#' @rdname sub_vector-methods
setMethod(f="sub_vector", signature=signature(x="FileBackedNeuroVec", i="numeric"),
          def=function(x, i) {
            mat <- read_mapped_vols(x@meta, i)
            sp <- add_dim(drop_dim(space(x)), length(i))
            DenseNeuroVec(mat, sp)
          })


#' as.list
#'
#' convert FileBackedNeuroVec to list of \code{\linkS4class{DenseNeuroVol}} objects.
#' @rdname as.list-methods
#' @export
setMethod(f="as.list", signature=signature(x = "FileBackedNeuroVec"), def=function(x) {
  D4 <- dim(x)[4]

  f <- function(i) {
    drop(sub_vector(x, i))
  }

  #deferred_list(lapply(seq(1, D4),
  #       function(i) { f } ))

  deflist::deflist(f, D4)

})

#' @noRd
setMethod(f="linear_access", signature=signature(x = "FileBackedNeuroVec", i = "numeric"),
          def=function (x, i) {
            read_mapped_data(x@meta, i)
          })


#' @noRd
setMethod(f="linear_access", signature=signature(x = "H5NeuroVec", i = "numeric"),
          def=function (x, i) {
            els <- x@obj[["data/elements"]]
            g <- index_to_grid(space(x), i)
            apply(g, 1, function(r) {
              els[r[1], r[2], r[3], r[4]]
            })
          })



#' @rdname as.matrix-methods
#' @export
setMethod(f="as.matrix", signature=signature(x = "FileBackedNeuroVec"), def=function(x) {
  as(x, "matrix")
})

#' @export
setAs(from="FileBackedNeuroVec", to="matrix",
      function(from) {
        len <- prod(from@meta@dims[1:3])
        t(series(from, seq(1, len)))
      })


#' @export
#' @rdname concat-methods
# setMethod(f="concat", signature=signature(x="FileBackedNeuroVec", y="FileBackedNeuroVec"),
#           def=function(x,y,...) {
#             if (!all(dim(x)[1:3] == dim(y)[1:3])) {
#               stop("cannot concatenate arguments with different spatial dimensions")
#             }
#             if (!all(spacing(x) == spacing(y))) {
#               stop("cannot concatenate arguments with different voxel spacing")
#             }
#
#             ndat <- rbind(x@data, y@data)
#             d1 <- dim(x)
#             d2 <- dim(y)
#
#             rest <- list(...)
#
#
#             if (length(rest) >= 1) {
#               mat <- do.call(rbind, map(rest, ~ .@data))
#
#               ndim <- c(d1[1:3], d1[4] + d2[4] + nrow(mat))
#               ndat <- rbind(ndat, mat)
#               nspace <- NeuroSpace(ndim, spacing(x@space),  origin(x@space), axes(x@space), trans(x@space))
#               SparseNeuroVec(ndat, nspace, mask=x@mask)
#             } else {
#               ndim <- c(d1[1:3], d1[4] + d2[4])
#               nspace <- NeuroSpace(ndim, spacing(x@space),  origin(x@space), axes(x@space), trans(x@space))
#               SparseNeuroVec(ndat, nspace, mask=x@mask)
#             }
#
#           })
