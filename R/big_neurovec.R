#' @keywords internal
BigNeuroVec <- function(data, space, mask,type = c("double", "float", "integer"), backingfile=tempfile()) {
  type <- match.arg(type)
  stopifnot(inherits(space, "NeuroSpace"))

  p <- prep_sparsenvec(data, space, mask)

  fbm <- bigstatsr::as_FBM(p$data, type=type, backingfile=backingfile)

  new("BigNeuroVec", space=p$space, mask=p$mask,
      map=IndexLookupVol(space(p$mask), as.integer(which(p$mask))), data=fbm)

}
