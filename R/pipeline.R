
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
  #deferred_list(replicate(len, f))
  ret <- structure(vector(mode="list", length=0), f=f, len=len, class=c("deferred_list2", "pairlist"))
  class(ret) <- c("deferred_list2", "list")
  ret
}


#' @export
#' @method print deferred_list2
print.deferred_list2 <- function(x,...) {
  cat("deferred_list2: ", attr(x, "len"), " elements. \n")
}

#' @export
#' @method as.list deferred_list2
as.list.deferred_list2 <- function(x,...) {
  purrr::map(seq(1,attr(x, "len")), ~ x[[.]])
}


#' @keywords internal
#' @export
`[[.deferred_list2` <- function (x, i)  {
  #ff <- NextMethod()
  #ff(i)
  #stopifnot(i <= x$len)
  stopifnot(i <= attr(x, "len") && i > 0)
  #x$f(i)
  attr(x, "f")(i)
}

#' @keywords internal
#' @export
`[.deferred_list2` <- function (x, i)  {
  #ff <- NextMethod()
  #lapply(seq_along(i), function(j) x$f(i[j]))
  f <- attr(x, "f")
  lapply(seq_along(i), function(j) f(i[j]))
}

#' @keywords internal
#' @export
length.deferred_list2 <- function (x)  {
  #x$len
  attr(x, "len")
}




#' #' @keywords internal
#' eval.deferred_list <- function(x) {
#'   lapply(1:length(x), function(i)
#'     x[[i]])
#' }
