#' @include all_class.R
{}
#' @include all_generic.R
{}



#' @export
#' @rdname ndim-methods
setMethod(f="ndim", signature=signature(x = "NeuroObj"),
          def=function(x) length(dim(x@space)))


#' dim of \code{NeuroObj} object
#' @param x the object
#' @rdname dim-methods
#' @export
setMethod(f="dim", signature=signature(x = "NeuroObj"),
          def=function(x) dim(x@space))


#' @export
#' @rdname space-methods
setMethod(f="space", signature=signature(x = "ROICoords"),
          def=function(x) x@space)

#' @export
#' @rdname spacing-methods
setMethod(f="spacing", signature=signature(x = "ROICoords"),
          def=function(x) spacing(x@space))

#' @export
#' @rdname space-methods
setMethod(f="space", signature=signature(x = "NeuroObj"),
          def=function(x) x@space)

#' @export
#' @rdname trans-methods
setMethod(f="trans", signature=signature(x = "NeuroObj"),
          def=function(x) trans(x@space))


#' @export
#' @rdname spacing-methods
setMethod(f="spacing",signature= signature(x = "NeuroObj"),
          def=function(x) {
            sp <- space(x)
            spacing(sp)
          })
