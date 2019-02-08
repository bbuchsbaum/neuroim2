
#' deferred_list
#' @keywords internal
#' @export
deferred_list <- function(fs) {
  assert_that(all(map_lgl(fs, is.function)))
  ret <- as.list(fs)
  structure(ret, class=c("deferred_list", "list"))
}


#' @export
#' @method print deferred_list
print.deferred_list <- function(x,...) {
  cat("deferred_list: ", length(x), " elements. \n")
}

#' @export
#' @method as.list deferred_list
as.list.deferred_list <- function(x,...) {
  map(seq_along(x), function(i) x[[i]])
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
  lapply(seq_along(i), function(j) ff[[j]](i[j]))

}




#' deferred_list2
#' @keywords internal
#' @export
deferred_list2 <- function(f, len=1) {
  assert_that(is.function(f))
  structure(f, len=len, class=c("deferred_list2"))
}


#' @export
#' @method print deferred_list
print.deferred_list2 <- function(x,...) {
  cat("deferred_list2: ", attr(x, "len"), " elements. \n")
}

#' @export
#' @method as.list deferred_list
as.list.deferred_list2 <- function(x,...) {
  purrr::map(seq(1,attr(x, "len")), ~ x[[.]])
}


#' @keywords internal
#' @export
`[[.deferred_list2` <- function (x, i)  {
  #ff <- NextMethod()
  #ff(i)
  stopifnot(i <= attr(x, "len"))
  x(i)
}

#' @keywords internal
#' @export
`[.deferred_list2` <- function (x, i)  {
  #ff <- NextMethod()
  lapply(seq_along(i), function(j) x(i[j]))
}



#' #' @keywords internal
#' eval.deferred_list <- function(x) {
#'   lapply(1:length(x), function(i)
#'     x[[i]])
#' }
