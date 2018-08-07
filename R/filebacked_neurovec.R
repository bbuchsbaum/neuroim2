
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


