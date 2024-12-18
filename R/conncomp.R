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
#' # Create a simple 3D binary mask with two components
#' mask <- array(FALSE, c(4, 4, 4))
#' mask[1:2, 1:2, 1:2] <- TRUE  # First component
#' mask[3:4, 3:4, 3:4] <- TRUE  # Second component
#'
#' # Extract components using 6-connectivity
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
#' @references
#' Rosenfeld, A., & Pfaltz, J. L. (1966). Sequential operations in digital 
#' picture processing. Journal of the ACM, 13(4), 471-494.
#'
#' @seealso 
#' \code{\link{array}} for creating 3D arrays
#' 
#' @importFrom purrr map flatten_dbl
#' @importFrom stats setNames
#'
#' @export
conn_comp_3D <- function(mask, connect = c("26-connect", "18-connect", "6-connect")) {
  # Input validation with more informative messages
  if (!is.array(mask) || length(dim(mask)) != 3) {
    stop("'mask' must be a 3D array")
  }
  if (!is.logical(mask[1])) {
    stop("'mask' must be a logical array")
  }

  connect <- match.arg(connect)

  # Use integer arrays for better memory efficiency
  nodes <- seq_len(length(mask) %/% 9)  # Initialize nodes with identity mapping
  labels <- array(0L, dim(mask))        # Use integer array
  
  DIM <- dim(mask)
  
  # Pre-compute neighborhood patterns more efficiently
  local.mask <- switch(connect,
    "6-connect" = {
      matrix(c(
        -1,0,0,  1,0,0,  0,-1,0,  0,1,0,  0,0,-1,  0,0,1
      ), ncol=3, byrow=TRUE)
    },
    "18-connect" = {
      matrix(c(
        -1,0,0,  1,0,0,   0,-1,0,  0,1,0,   0,0,-1,  0,0,1,  # 6-neighbors
        -1,-1,0, -1,1,0,  1,-1,0,  1,1,0,                     # edge neighbors xy
        -1,0,-1, -1,0,1,  1,0,-1,  1,0,1,                     # edge neighbors xz
        0,-1,-1, 0,-1,1,  0,1,-1,  0,1,1                      # edge neighbors yz
      ), ncol=3, byrow=TRUE)
    },
    # 26-connect: more efficient than expand.grid
    matrix(c(
      -1,-1,-1, -1,-1,0, -1,-1,1,  -1,0,-1, -1,0,0, -1,0,1,  -1,1,-1, -1,1,0, -1,1,1,
       0,-1,-1,  0,-1,0,  0,-1,1,   0,0,-1,         0,0,1,    0,1,-1,  0,1,0,  0,1,1,
       1,-1,-1,  1,-1,0,  1,-1,1,   1,0,-1,  1,0,0, 1,0,1,    1,1,-1,  1,1,0,  1,1,1
    ), ncol=3, byrow=TRUE)
  )
  
  # Transpose once for efficiency
  tlocal.mask <- t(local.mask)
  
  # Optimized neighbor function with fewer allocations
  neighbors <- function(vox) {
    # Add current voxel to each neighbor offset
    vox.hood <- t(tlocal.mask + vox)
    
    # Handle boundary cases more efficiently
    if (any(vox == 1) || any(vox == DIM)) {
      valid <- vox.hood[,1] >= 1 & vox.hood[,1] <= DIM[1] &
               vox.hood[,2] >= 1 & vox.hood[,2] <= DIM[2] &
               vox.hood[,3] >= 1 & vox.hood[,3] <= DIM[3]
      vox.hood <- vox.hood[valid, , drop=FALSE]
    }
    
    # Get labeled neighbors
    has_label <- labels[cbind(vox.hood[,1], vox.hood[,2], vox.hood[,3])] != 0
    vox.hood[has_label, , drop=FALSE]
  }
  
  # Path compression in find() for better performance
  find <- function(i) {
    root <- i
    # Find root
    while (nodes[root] != root) {
      root <- nodes[root]
    }
    # Path compression
    while (nodes[i] != root) {
      parent <- nodes[i]
      nodes[i] <- root
      i <- parent
    }
    root
  }
  
  # First pass: Initial labeling
  nextlabel <- 1L
  grid <- which(mask > 0, arr.ind=TRUE)  # More efficient than .indexToGrid
  
  for (i in seq_len(nrow(grid))) {
    vox <- grid[i,]
    nabes <- neighbors(vox)
    
    if (nrow(nabes) == 0) {
      nodes[nextlabel] <- nextlabel
      labels[vox[1], vox[2], vox[3]] <- nextlabel
    } else {
      L <- labels[cbind(nabes[,1], nabes[,2], nabes[,3])]
      ML <- min(L)
      labels[vox[1], vox[2], vox[3]] <- ML
      nodes[nextlabel] <- ML
      
      # Union operation with unique labels for efficiency
      for (lab in unique(L)) {
        rootx <- find(lab)
        nodes[rootx] <- find(ML)
      }
    }
    
    nextlabel <- nextlabel + 1L
  }
  
  # Second pass: Resolve labels more efficiently
  label_indices <- which(labels > 0)
  labels[label_indices] <- vapply(labels[label_indices], find, integer(1))
  
  # Create output volumes efficiently
  labs <- labels[label_indices]
  clusters <- sort(table(labs), decreasing=TRUE)
  
  # Handle case of no components
  if (length(clusters) == 0) {
    return(list(
      index = array(0L, DIM),
      size = array(0L, DIM)
    ))
  }
  
  # Size volume
  SVol <- array(0L, DIM)
  SVol[label_indices] <- clusters[as.character(labs)]
  
  # Index volume with consecutive indices
  indices <- seq_along(clusters)
  names(indices) <- names(clusters)
  IVol <- array(0L, DIM)
  IVol[label_indices] <- indices[as.character(labs)]
  
  list(index=IVol, size=SVol)
}
