
#' deferred_list
#' @keywords internal
#' @export
deferred_list <- function(fs) {
  stopifnot(all(map_lgl(fs, is.function)))
  ret <- as.list(fs)
  structure(ret, class="deferred_list")
}



#' @keywords internal
#' @export
`[[.deferred_list` <- function (x, i)  {
 ff <- NextMethod()

 ff(i)
}

#' @keywords internal
#' @export
`[.deferred_list` <- function (x, i)  {
  ff <- NextMethod()
  deferred_list(ff)
}



#' #' @keywords internal
#' eval.deferred_list <- function(x) {
#'   lapply(1:length(x), function(i)
#'     x[[i]])
#' }
