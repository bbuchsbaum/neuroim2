#' @include all_class.R
{}
#' @include all_generic.R
{}

#' Two-Dimensional Neuroimaging Data Slice
#'
#' @title NeuroSlice: 2D Neuroimaging Data Container
#'
#' @description
#' Creates a \code{NeuroSlice} object representing a two-dimensional slice of neuroimaging data
#' with associated spatial information. This class is particularly useful for working with
#' individual slices from volumetric neuroimaging data or for visualizing 2D cross-sections.
#'
#' @param data A vector or matrix containing the slice data values.
#' @param space An object of class \code{\linkS4class{NeuroSpace}} defining the spatial
#'   properties (dimensions, spacing, origin) of the slice.
#' @param indices Optional integer vector. When \code{data} is provided as a 1D vector,
#'   \code{indices} specifies the linear indices where the data values should be placed
#'   in the 2D slice. Useful for creating sparse slices. Default is NULL.
#'
#' @return A new object of class \code{\linkS4class{NeuroSlice}}.
#'
#' @section Input Validation:
#' The function performs several validation checks:
#' \itemize{
#'   \item Verifies that \code{space} is 2-dimensional
#'   \item Ensures data dimensions are compatible with \code{space}
#'   \item Validates \code{indices} when provided for sparse initialization
#' }
#'
#' @section Data Handling:
#' The function supports two initialization modes:
#' \itemize{
#'   \item Dense mode (indices = NULL):
#'     \itemize{
#'       \item Data is reshaped if necessary to match space dimensions
#'       \item Dimensions must match exactly after reshaping
#'     }
#'   \item Sparse mode (indices provided):
#'     \itemize{
#'       \item Creates a zero-initialized matrix matching space dimensions
#'       \item Places data values at specified indices
#'     }
#' }
#'
#' @examples
#' # Create a 64x64 slice space
#' slice_space <- NeuroSpace(c(64, 64), spacing = c(2, 2))
#'
#' # Example 1: Dense slice from matrix
#' slice_data <- matrix(rnorm(64*64), 64, 64)
#' dense_slice <- NeuroSlice(slice_data, slice_space)
#'
#' # Example 2: Dense slice from vector
#' vec_data <- rnorm(64*64)
#' vec_slice <- NeuroSlice(vec_data, slice_space)
#'
#' # Example 3: Sparse slice with specific values
#' n_points <- 100
#' sparse_data <- rnorm(n_points)
#' sparse_indices <- sample(1:(64*64), n_points)
#' sparse_slice <- NeuroSlice(sparse_data, slice_space, indices = sparse_indices)
#'
#' @seealso
#' \code{\linkS4class{NeuroSpace}} for defining spatial properties,
#' \code{\linkS4class{NeuroVol}} for 3D volumetric data,
#' \code{\link{plot}} for visualization methods
#'
#' @export
NeuroSlice <- function(data, space, indices = NULL) {
  assert_that(ndim(space) == 2,
              msg = "Space must be 2-dimensional for NeuroSlice")

  if (is.null(indices)) {
    if (length(dim(data)) != 2) {
      assert_that(length(data) == prod(dim(space)[1:2]),
                  msg = "Data length must match space dimensions")
      data <- matrix(data, dim(space)[1], dim(space)[2])
    }

    assert_that(all(dim(data) == dim(space)),
                msg = "Data dimensions must match space dimensions")
    new("NeuroSlice", .Data=data, space=space)

  } else {
    mdat <- matrix(0, dim(space))
    mdat[indices] <- data
    new("NeuroSlice", .Data=mdat, space=space)
  }
}


#' Convert Grid Coordinates to Linear Indices
#'
#' @title Grid to Linear Index Conversion
#' @description
#' Converts 2D grid coordinates to linear indices for a \code{NeuroSlice} object.
#'
#' @param x A \code{NeuroSlice} object
#' @param coords Either a numeric vector of length 2 or a matrix with 2 columns,
#'   representing (x,y) coordinates in the slice grid
#'
#'
#' @rdname grid_to_index-methods
#'
#' @examples
#' slice_space <- NeuroSpace(c(10, 10))
#' slice_data <- matrix(1:100, 10, 10)
#' slice <- NeuroSlice(slice_data, slice_space)
#'
#' # Convert single coordinate
#' idx <- grid_to_index(slice, c(5, 5))
#'
#' # Convert multiple coordinates
#' coords <- matrix(c(1,1, 2,2, 3,3), ncol=2, byrow=TRUE)
#' indices <- grid_to_index(slice, coords)
#'
#' @seealso \code{\link{index_to_grid}} for the inverse operation
#'
#' @export
setMethod(f="grid_to_index",
          signature=signature(x = "NeuroSlice", coords="matrix"),
          def=function(x, coords) {
            callGeneric(x@space, coords)
          })

#' @rdname grid_to_index-methods
#' @export
setMethod(f="grid_to_index",
          signature=signature(x = "NeuroSlice", coords="numeric"),
          def=function(x, coords) {
            callGeneric(x@space, coords)
          })

#' Convert Linear Indices to Grid Coordinates
#'
#' @title Linear Index to Grid Coordinate Conversion
#' @description
#' Converts linear indices to 2D grid coordinates for a \code{NeuroSlice} object.
#'
#' @param x A \code{NeuroSlice} object
#' @param idx Integer vector of linear indices to convert
#'
#'
#' @examples
#' slice_space <- NeuroSpace(c(10, 10))
#' slice_data <- matrix(1:100, 10, 10)
#' slice <- NeuroSlice(slice_data, slice_space)
#'
#' # Convert single index
#' coords <- index_to_grid(slice, 55)
#'
#' # Convert multiple indices
#' indices <- c(1, 25, 50, 75, 100)
#' coords_mat <- index_to_grid(slice, indices)
#'
#' @seealso \code{\link{grid_to_index}} for the inverse operation
#'
#' @export
#' @rdname index_to_grid-methods
setMethod(f="index_to_grid",
          signature=signature(x = "NeuroSlice", idx="numeric"),
          def=function(x, idx) {
            callGeneric(x@space, idx)
          })

#' Plot a NeuroSlice
#'
#' @name plot,NeuroSlice-method
#' @param cmap Color map to use for plotting, defaults to grayscale
#' @param irange Intensity range for scaling the plot values, defaults to the data range
#' @param legend Logical indicating whether to display the color legend. Defaults to TRUE.
#'
#' @return a ggplot2 object
#'
#' @details
#' The plot method uses \code{ggplot2} to create a raster visualization of the slice data.
#' The intensity values are mapped to colors using the specified colormap and range.
#'
#' @details when `x` is a NeuroSlice object, the plot method returns a \code{ggplot2} object containing the raster visualization of the slice data.
#'         The plot can be further customized using standard ggplot2 functions.
#'
#' @examples
#' # Create example slice
#' slice_space <- NeuroSpace(c(100, 100))
#' slice_data <- matrix(rnorm(100*100), 100, 100)
#' slice <- NeuroSlice(slice_data, slice_space)
#' \donttest{
#' # Basic plot
#' plot(slice)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn xlab ylab theme_bw
#' @importFrom grDevices gray
#' @export
#' @rdname plot-methods
setMethod("plot",
          signature=signature(x="NeuroSlice"),
          def=function(x, cmap=gray(seq(0,1,length.out=255)),
                      irange=range(x, na.rm=TRUE),
                      legend=TRUE) {
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            {y = value = NULL}

            cds <- index_to_coord(space(x), 1:length(x))
            df1 <- data.frame(x = cds[,1], y = cds[,2], value = as.vector(x))

            ggplot2::ggplot(df1, ggplot2::aes(x = x, y = y, fill = value)) +
              ggplot2::geom_raster() +
              ggplot2::scale_fill_gradientn(colours = cmap,
                                           limits = irange,
                                           guide = if (legend) "colourbar" else "none") +
              ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::theme_bw()

          })


#' @rdname show-methods
#' @export
setMethod(f="show",
          signature=signature("NeuroSlice"),
          def=function(object) {
            # Get space information
            sp <- space(object)

            # Calculate statistics
            val_range <- range(object, na.rm=TRUE)
            n_na <- sum(is.na(object))
            mem_size <- format(object.size(object), units="auto")

            # Header
            cat("\n")
            cat(bold(blue("=== NeuroSlice Object ===")), "\n\n")

            # Type and Dimensions
            cat(bold(yellow("* Basic Information")), "\n")
            cat("  ", silver("Type:"), " ", class(object), "\n", sep="")
            cat("  ", silver("Dimensions:"), " ",
                paste(dim(object), collapse=" x "),
                " (", green(mem_size), ")", "\n", sep="")

            # Value Range and Stats
            cat("\n", bold(yellow("* Data Properties")), "\n", sep="")
            cat("  ", silver("Value Range:"), " [",
                blue(sprintf("%.2f", val_range[1])), ", ",
                blue(sprintf("%.2f", val_range[2])), "]", "\n", sep="")
            if (n_na > 0) {
                cat("  ", silver("Missing Values:"), " ",
                    red(n_na), " (",
                    sprintf("%.1f%%", 100*n_na/length(object)), ")",
                    "\n", sep="")
            }

            # Spatial Properties
            cat("\n", bold(yellow("* Spatial Properties")), "\n", sep="")
            cat("  ", silver("Spacing:"), " ",
                paste(sprintf("%.2f", sp@spacing), collapse=" x "),
                "\n", sep="")
            cat("  ", silver("Origin:"), "  ",
                paste(sprintf("%.2f", sp@origin), collapse=" x "),
                "\n", sep="")
            cat("  ", silver("Axes:"), "    ",
                green(sp@axes@i@axis), " x ",
                green(sp@axes@j@axis), "\n", sep="")

            # Footer
            cat("\n", blue("=" = 28), "\n", sep="")
          })

#' @import assertthat
#' @keywords internal
#' @importFrom grDevices col2rgb gray heat.colors
#' @noRd
mapToColors <- function (imslice, col = heat.colors(128, alpha = 1), zero_col = "#00000000",
                         alpha = 1, irange = range(imslice), threshold = c(0, 0)) {

  assertthat::assert_that(diff(irange) >= 0)
  drange <- diff(irange)
  mcols <- (imslice - irange[1])/diff(irange) * (length(col) -1) + 1
  mcols[mcols < 1] <- 1
  mcols[mcols > length(col)] <- length(col)
  imcols <- col[mcols]

  if (!is.vector(imslice)) {
    dim(imcols) <- dim(imslice)
  }

  imcols[imslice == 0] <- zero_col

  if (diff(threshold) > 0) {
    imcols[(imslice >= threshold[1]) & (imslice <= threshold[2])] <- "#00000000"
  }

  if (alpha < 1) {
    rgbmat <- col2rgb(imcols, alpha = TRUE)
    rgbmat <- rgbmat/255
    rgbmat[4, ] <- rgbmat[4, ] * alpha

    if (is.vector(imslice)) {
      array(t(rgbmat), c(length(imslice), 4))
    } else {
      array(t(rgbmat), c(dim(imslice), 4))
    }
  }
  else {
    imcols
  }
}

#' @rdname mask-methods
#' @export
setMethod("mask", "NeuroSlice",
          function(x) {
            # Create a 3D mask with the slice dimension expanded to 1
            dims <- dim(x)
            
            # NeuroSlice is always 2D, expand to 3D for mask
            if (length(dims) == 2) {
              # Add a third dimension of size 1
              mask_dims <- c(dims[1], dims[2], 1)
              # Use default spacing and origin since NeuroSlice doesn't have these
              LogicalNeuroVol(array(TRUE, mask_dims),
                            NeuroSpace(mask_dims))
            } else {
              # Shouldn't happen for NeuroSlice, but handle gracefully
              LogicalNeuroVol(array(TRUE, dims), 
                            NeuroSpace(dims))
            }
          })
