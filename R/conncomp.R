#' Extract Connected Components from a 3D Binary Mask
#'
#' @description
#' Identifies and labels connected components in a 3D binary mask using a two-pass algorithm.
#' The function supports different connectivity constraints and returns both component indices
#' and their sizes.
#'
#' @param mask A 3D logical array representing the binary mask
#' @param connect A character string specifying the connectivity constraint. One of 
#'   "26-connect" (default), "18-connect", or "6-connect"
#'
#' @return A list with the following components:
#'   \item{index}{A 3D array of integers. Each non-zero value represents the cluster index 
#'                of the connected component for that voxel. Zero values indicate background.}
#'   \item{size}{A 3D array of integers. Each non-zero value represents the size (number of voxels) 
#'               of the connected component that the voxel belongs to. Zero values indicate background.}
#'
#' @details
#' The function implements an efficient two-pass connected component labeling algorithm:
#'
#' \itemize{
#'   \item First pass: Assigns provisional labels and builds an equivalence table using
#'         a union-find data structure for label resolution
#'   \item Second pass: Resolves label conflicts and assigns final component labels
#' }
#'
#' The connectivity options determine which voxels are considered adjacent:
#' \itemize{
#'   \item 6-connect: Only face-adjacent voxels (±1 step along each axis)
#'   \item 18-connect: Face and edge-adjacent voxels
#'   \item 26-connect: Face, edge, and vertex-adjacent voxels (all neighbors in a 3x3x3 cube)
#' }
#'
#' Time complexity is O(n) where n is the number of voxels in the mask, with
#' additional O(k) space for the union-find data structure where k is the number
#' of provisional labels.
#'
#' @examples
#' # Create a simple 3D binary mask with two disconnected components
#' mask <- array(FALSE, c(4, 4, 4))
#' mask[1:2, 1:2, 1:2] <- TRUE  # First component
#' mask[3:4, 3:4, 3:4] <- TRUE  # Second component
#'
#' # Extract components using different connectivity patterns
#' comps <- conn_comp_3D(mask, connect = "6-connect")
#'
#' # Number of components
#' max_comps <- max(comps$index)
#' cat("Found", max_comps, "components\n")
#'
#' # Size of each component
#' unique_sizes <- unique(comps$size[comps$size > 0])
#' cat("Component sizes:", paste(unique_sizes, collapse=", "), "\n")
#'
#' # Try with different connectivity
#' comps_26 <- conn_comp_3D(mask, connect = "26-connect")
#' cat("Number of components with 26-connectivity:", max(comps_26$index), "\n")
#'
#' @references
#' Rosenfeld, A., & Pfaltz, J. L. (1966). Sequential operations in digital 
#' picture processing. Journal of the ACM, 13(4), 471-494.
#'
#' @seealso 
#' \code{\link[base]{array}} for creating 3D arrays,
#' \code{\link{ClusteredNeuroVol}} for working with clustered neuroimaging data
#' 
#' @importFrom purrr map flatten_dbl
#' @importFrom stats setNames
#'
#' @export
conn_comp_3D <- function(mask, connect = c("26-connect", "18-connect", "6-connect")) {
  # Input validation with more informative messages
  if (!is.array(mask) || length(dim(mask)) != 3) {
    cli::cli_abort("{.arg mask} must be a 3D array.")
  }
  if (!is.logical(mask)) {
    cli::cli_abort("{.arg mask} must be a logical array.")
  }

  connect <- match.arg(connect)
  connectivity <- switch(connect, "6-connect" = 6L, "18-connect" = 18L, "26-connect" = 26L)

  DIM <- dim(mask)
  # The labelling, the union-find and the size ordering all happen in one
  # compiled pass. The interpreted version this replaces was ~1000x slower and
  # got worse with volume size; see src/conncomp.cpp for why the numbering is
  # reproduced exactly rather than simplified.
  res <- conn_comp_labels_cpp(as.vector(mask), as.integer(DIM), connectivity)

  list(index = array(res$index, DIM), size = array(res$size, DIM))
}
