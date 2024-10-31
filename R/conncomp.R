#' Extract Connected Components from a 3D Mask
#'
#' @description
#' This function identifies and labels connected components in a 3D binary mask using
#' different connectivity constraints.
#'
#' @param mask A 3D binary array representing the mask. Must be a logical array.
#' @param connect A character string specifying the connectivity constraint. 
#'   Options are "6-connect", "18-connect", or "26-connect". Default is "26-connect".
#'
#' @return A list with two elements:
#'   \item{index}{A 3D array containing the cluster index of the connected component for each voxel.}
#'   \item{size}{A 3D array containing the size of the connected component for each voxel.}
#'
#' @details
#' The function implements a two-pass algorithm to identify connected components:
#' \itemize{
#'   \item First pass: Assigns initial labels and establishes connectivity.
#'   \item Second pass: Resolves label conflicts and assigns final labels.
#' }
#' 
#' The connectivity options determine which voxels are considered neighbors:
#' \itemize{
#'   \item 6-connect: Only face-adjacent voxels
#'   \item 18-connect: Face and edge-adjacent voxels
#'   \item 26-connect: Face, edge, and vertex-adjacent voxels
#' }
#'
#' @examples
#' # Generate a random 3D binary mask
#' set.seed(123)
#' dat <- array(as.logical(runif(10*10*10) > 0.5), c(10, 10, 10))
#' 
#' # Extract connected components using different connectivity constraints
#' res1 <- conn_comp_3D(dat, connect = "6-connect")
#' res2 <- conn_comp_3D(dat, connect = "18-connect")
#' res3 <- conn_comp_3D(dat, connect = "26-connect")
#' 
#' # Compare the number of connected components
#' cat("Number of components (6-connect):", max(res1$index), "\n")
#' cat("Number of components (18-connect):", max(res2$index), "\n")
#' cat("Number of components (26-connect):", max(res3$index), "\n")
#'
#' @seealso 
#' \code{\link{array}} for creating 3D arrays
#' 
#' @importFrom purrr map flatten_dbl
#'
#' @export
conn_comp_3D <- function(mask, connect = c("26-connect", "18-connect", "6-connect")) {
	stopifnot(length(dim(mask)) == 3 && is.logical(mask[1]))

  connect <- match.arg(connect)

	nodes <- numeric(length(mask)/9)
	labels <- array(0, dim(mask))

	DIM <- dim(mask)
	xdim <- DIM[1]
	ydim <- DIM[2]
	zdim <- DIM[3]

	local.mask <- if (connect == "6-connect") {
	  as.matrix(
	    rbind(expand.grid(x=c(-1,0,1), y=0, z=0),
	          expand.grid(x=0, y=c(-1, 1), z=0),
	          expand.grid(x=0, y=0, z=c(-1, 1)))
	  )
	} else if (connect == "18-connect") {
	  as.matrix(rbind(
	    expand.grid(x=c(-1,0,1), y=0, z=0),
	    expand.grid(x=0, y=c(-1,1), z=0),
	    expand.grid(x=0, y=0, z=c(-1,1)),
	    expand.grid(x=c(-1,1), y=c(-1,1), z=0),
	    expand.grid(x=c(-1,1), y=0, z=c(-1,1)),
	    expand.grid(x=0, y=c(-1,1), z=c(-1,1)))
	  )
	} else {
	  as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1), z=c(-1,0,1)))
	}

	dimnames(local.mask) <- NULL
	local.mask <- local.mask[-(ceiling(nrow(local.mask)/2)),,drop=FALSE]
	tlocal.mask <- t(local.mask)

	neighbors <- function(vox) {

	  vox.hood <- t(tlocal.mask + vox)
	  if (any(vox == 1) || any(vox == DIM)) {
	    vox.hood <- vox.hood[apply(vox.hood, 1, function(coords) {
	      all(coords > 1 & coords <= DIM)
	    }),,drop=FALSE]
	  }

	  vox.hood[labels[vox.hood] != 0,,drop=F]
	}


	find <- function(i) {
		while (nodes[i] != i) {
			i <- nodes[i]
		}

		nodes[i]
	}


	nextlabel <- 1

	grid <-  .indexToGrid(which(mask>0), dim(mask))

	for (i in 1:NROW(grid)) {
		vox <- grid[i,]
		nabes <- neighbors(vox)
		if (nrow(nabes) == 0) {
			nodes[nextlabel] <- nextlabel
			labels[vox[1],vox[2],vox[3]] <- nextlabel
		} else {
			L <- labels[nabes]
			ML <- min(L)
			labels[vox[1],vox[2], vox[3]] <- ML
			nodes[nextlabel] <- ML
			for (lab in L) {
				rootx <- find(lab)
				nodes[rootx] <- find(ML)
			}
		}

		nextlabel <- nextlabel + 1
	}



	## pass2
	for (k in 1:zdim) {
		for (j in 1:ydim) {
			for (i in 1:xdim) {
				if (labels[i,j,k] > 0) {
					labels[i,j,k] <- find(labels[i,j,k])
				}
			}
		}
	}

	labs <- labels[labels!=0]
	forelabs <- labels > 0

	clusters <- sort(table(labs), decreasing=TRUE)
	SVol <- array(0, dim(mask))
	SVol[forelabs] <- clusters[as.character(labs)]

	indices <- 1:length(clusters)
	names(indices) <- names(clusters)
	IVol <- array(0, dim(mask))
	IVol[forelabs] <- indices[as.character(labs)]


	list(index=IVol, size=SVol)

}




