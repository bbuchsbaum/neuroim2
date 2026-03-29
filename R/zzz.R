#' @useDynLib neuroim2, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats sd
NULL

.onUnload <- function(libpath) {
  if (exists(".clear_mmap_cache", mode = "function")) {
    .clear_mmap_cache()
  }
  invisible(NULL)
}
