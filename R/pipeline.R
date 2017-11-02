#' @keywords internal
deferred_list <- function(fs) {
  stopifnot(all(sapply(fs, is.function)))
  ret <- as.list(fs)
  structure(ret, class="deferred_list")
}



#' @keywords internal
`[[.deferred_list` <- function (x, i)  {
 ff <- NextMethod()
 ff(i)
}

#' #' @keywords internal
#' eval.deferred_list <- function(x) {
#'   lapply(1:length(x), function(i)
#'     x[[i]])
#' }
