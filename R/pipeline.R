
#' deferred_list
#' @keywords internal
#' @export
deferred_list <- function(fs) {
  assert_that(all(map_lgl(fs, is.function)))
  ret <- as.list(fs)
  structure(ret, class=c("deferred_list", "list"))
}


#' @export
print.deferred_list <- function(x) {
  cat("deferred_list: ", length(x), " elements. \n")
}

#' @export
as.list.deferred_list <- function(x) {
  map(1:length(x), function(i) x[[i]])
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
