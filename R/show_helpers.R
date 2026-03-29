#' @include all_class.R
NULL

#' @importFrom cli format_inline rule
#' @importFrom utils object.size
NULL

#' Internal helpers for consistent show() formatting
#'
#' All show() methods in the package should use these helpers to ensure
#' a uniform visual style. The output uses \pkg{cli} for section rules
#' and plain \code{cat()} for field lines.
#'
#' @name show-helpers
#' @keywords internal
NULL

#' Print a show() header line
#'
#' @param class_name Character class name.
#' @param tag Optional tag shown in brackets after the class name.
#' @keywords internal
#' @noRd
show_header <- function(class_name, tag = NULL) {
  if (!is.null(tag) && nzchar(tag)) {
    cat(cli::format_inline("{.cls {class_name}} [{tag}]"), "\n")
  } else {
    cat(cli::format_inline("{.cls {class_name}}"), "\n")
  }
}

#' Print a section rule
#'
#' @param title Section title.
#' @keywords internal
#' @noRd
show_rule <- function(title) {
  cat(cli::rule(left = title), "\n")
}

#' Print a label : value field
#'
#' @param label Field name (left-padded to \code{width}).
#' @param ... Pasted together to form the value string.
#' @param width Padding width for the label column.
#' @keywords internal
#' @noRd
show_field <- function(label, ..., width = 14) {
  lbl <- formatC(label, width = -width, flag = "-")
  val <- paste0(...)
  cat("  ", lbl, ": ", val, "\n", sep = "")
}

#' Safely retrieve orientation axis codes
#'
#' @param sp A \code{NeuroSpace} object.
#' @return A character string like \code{"RAS"}, or \code{"N/A"} on failure.
#' @keywords internal
#' @noRd
safe_axcodes <- function(sp) {
  tryCatch(
    paste(affine_to_axcodes(trans(sp)), collapse = ""),
    error = function(e) "N/A"
  )
}

#' Format memory size as a human-readable string
#'
#' @param obj Any R object.
#' @return Character string like \code{"1.2 MB"}.
#' @keywords internal
#' @noRd
format_mem <- function(obj) {
  format(utils::object.size(obj), units = "auto")
}
