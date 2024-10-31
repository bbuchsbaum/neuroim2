#' Create a NeuroHyperVec instance
#'
#' @param ... One or more instances of type NeuroVec
#' @return A NeuroHyperVec object
#' @export
#'
#' @examples
#' space <- NeuroSpace(dim = c(64, 64, 32), origin = c(0, 0, 0), spacing = c(3, 3, 4))
#' vec1 <- NeuroVec(data = array(rnorm(64*64*32*100), dim = c(64, 64, 32, 100)), space = space)
#' vec2 <- NeuroVec(data = array(rnorm(64*64*32*100), dim = c(64, 64, 32, 100)), space = space)
#' hyper_vec <- NeuroHyperVec(vec1, vec2)
NeuroHyperVec <- function(...) {
  vecs <- list(...)
  assertthat::assert_that(all(purrr::map_lgl(vecs, ~ inherits(., "NeuroVec"))))
  
  # Check consistency of spatial dimensions and spacing
  first_vec <- vecs[[1]]
  assertthat::assert_that(all(purrr::map_lgl(vecs, ~ all(dim(.)[1:3] == dim(first_vec)[1:3]))))
  assertthat::assert_that(all(purrr::map_lgl(vecs, ~ all(spacing(.) == spacing(first_vec)))))
  assertthat::assert_that(all(purrr::map_lgl(vecs, ~ dim(.)[4] == dim(first_vec)[4])))
  
  # Create the NeuroHyperVec object
  new("NeuroHyperVec", space = space(first_vec), vecs = vecs)
}

#' Get dimensions of a NeuroHyperVec object
#'
#' @param x A NeuroHyperVec object
#' @return A vector of dimensions
#' @export
setMethod("dim", "NeuroHyperVec", function(x) {
  c(dim(x@vecs[[1]]), length(x@vecs))
})

#' Get the number of NeuroVec objects in a NeuroHyperVec
#'
#' @param x A NeuroHyperVec object
#' @return The number of NeuroVec objects
#' @export
setMethod("length", "NeuroHyperVec", function(x) {
  length(x@vecs)
})


#' Get a single NeuroVec from a NeuroHyperVec
#'
#' @param x A NeuroHyperVec object
#' @param i Index of the NeuroVec to retrieve
#' @return A NeuroVec object
#' @export
setMethod("[[", "NeuroHyperVec", function(x, i) {
  x@vecs[[i]]
})

#' Print method for NeuroHyperVec
#'
#' @param object A NeuroHyperVec object
#' @export
setMethod("show", "NeuroHyperVec", function(object) {
  cat("NeuroHyperVec object\n")
  cat("Dimensions:", paste(dim(object), collapse = " x "), "\n")
  cat("Number of NeuroVec objects:", length(object), "\n")
  cat("Spatial dimensions:", paste(dim(object@vecs[[1]])[1:3], collapse = " x "), "\n")
  cat("Time points per NeuroVec:", dim(object@vecs[[1]])[4], "\n")
})


