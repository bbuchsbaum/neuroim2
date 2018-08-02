#' CachedNeuroVecSource
#'
#' constructs a CachedNeuroVecSource object
#'
#' @param meta_info an object of class \code{\linkS4class{MetaInfo}}
#' @param indices a vector of 1D indices
#' @param mask a 3D \code{array} of type \code{logical}
#' @param bucket_vol a user supplied \code{\linkS4class{ClusteredNeuroVol}} defining the cache buckets
#' @rdname CachedSparseNeuroVecSource-class
CachedSparseNeuroVecSource <- function(meta_info, indices, bucket_vol) {

  stopifnot(length(dim(meta_info)) >= 3)
  stopifnot(all(indices >= 1 & indices <= dim(meta_info)[4]))

  D <- dim(meta_info)[1:3]
  mask <- bucket_vol@mask

  stopifnot(all(dim(mask) == D))

  new("CachedSparseNeuroVecSource", meta_info=meta_info, indices=indices, mask=mask, bucket_vol=bucket_vol)
}
