#'
#' neuroim2: neuroimaging data structures for analysis
#'
#' @description
#' The neuroim2 package provides tools and functions for analyzing and
#' manipulating neuroimaging data. It supports various neuroimaging formats
#' and offers a range of analysis techniques.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{read_vol}}: Read neuroimaging volumes
#'   \item \code{\link{write_vol}}: Write neuroimaging volumes
#'   \item \code{\link{NeuroVol}}: Create NeuroVol objects
#'   \item \code{\link{NeuroVec}}: Create NeuroVec objects
#' }
#'
#' @docType package
#' @name neuroim2-package
#' @aliases neuroim2
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom magrittr %>%
"_PACKAGE"

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
