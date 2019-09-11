#' Create an \code{NeuroHyperVec} instance from a variable length list of \code{NeuroVec} objects.
#'
#' @param ... one or more instances of type \code{NeuroVec}
#' @export
#'
#' @examples
#'
#' v1 <- NeuroVec(array(0,c(5,5,5,2)), space=NeuroSpace(dim=c(5,5,5,2)))
#' v2 <- NeuroVec(array(0,c(5,5,5,4)), space=NeuroSpace(dim=c(5,5,5,4)))
#' v3 <- NeuroVec(array(0,c(5,5,5,6)), space=NeuroSpace(dim=c(5,5,5,6)))
#' hv <- NeuroHyperVec(v1,v2,v3)
NeuroHyperVec <- function(...) {
  vecs <- list(...)
  assert_that(all(map_lgl(vecs, ~ inherits(., "NeuroVec"))))
  sp <- space(vecs[[1]])
  lens <- map_dbl(vecs, function(x) dim(x)[4])
  assert_that(length(lens) > 1)
  assert_that(all(lens == lens[1]))

  sp <- add_dim(sp, lens[1])

  new("NeuroHyperVec", space=sp, vecs=vecs)
}


#' @rdname series-methods
#' @export
setMethod("series", signature(x="NeuroHyperVec", i="ROICoords"),
          def=function(x,i) {
            callGeneric(x, coords(i))
          })

