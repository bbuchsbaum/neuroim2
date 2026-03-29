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
  if (ndim(space) != 2) {
    cli::cli_abort("Space must be 2-dimensional for {.cls NeuroSlice}, not {ndim(space)}D.")
  }

  if (is.null(indices)) {
    if (length(dim(data)) != 2) {
      if (length(data) != prod(dim(space)[1:2])) {
        cli::cli_abort("Data length ({length(data)}) must match space dimensions ({prod(dim(space)[1:2])}).")
      }
      data <- matrix(data, dim(space)[1], dim(space)[2])
    }

    if (!all(dim(data) == dim(space))) {
      cli::cli_abort("Data dimensions ({.val {dim(data)}}) must match space dimensions ({.val {dim(space)}}).")
    }
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
#' @rdname grid_to_index-methods
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
#' @aliases plot,NeuroSlice,ANY-method
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

            # Use pixel-grid coordinates so geom_raster stays memory-safe even
            # for oblique/sheared affine spaces where world coordinates are not
            # axis-aligned on a regular raster grid.
            df1 <- slice_df(x)

            ggplot2::ggplot(df1, ggplot2::aes(x = x, y = y, fill = value)) +
              ggplot2::geom_raster() +
              ggplot2::scale_fill_gradientn(colours = cmap,
                                           limits = irange,
                                           guide = if (legend) ggplot2::guide_colourbar(barheight = grid::unit(3, "cm")) else "none") +
              ggplot2::xlab("") + ggplot2::ylab("") +
              coord_neuro_fixed() +
              ggplot2::theme_bw()

          })


#' @rdname show-methods
#' @export
setMethod(f="show",
          signature=signature("NeuroSlice"),
          def=function(object) {
            sp <- space(object)
            show_header("NeuroSlice", format_mem(object))
            show_rule("Spatial")
            show_field("Dimensions", paste(dim(object), collapse = " x "))
            show_field("Spacing", paste(round(sp@spacing, 3), collapse = " x "))
            show_field("Origin", paste(round(sp@origin, 2), collapse = ", "))
            show_rule("Data")
            rng <- range(object, na.rm = TRUE)
            show_field("Range", sprintf("[%.3f, %.3f]", rng[1], rng[2]))
            n_na <- sum(is.na(object))
            if (n_na > 0) show_field("NAs", format(n_na, big.mark = ","))

            # Footer
            cat("\n", blue("=" = 28), "\n", sep="")
          })

#' Map intensity values to colors
#'
#' Convert intensity values (e.g., a 2D slice) into a color representation for
#' plotting and overlays.
#'
#' @param imslice A numeric vector or array of intensities.
#' @param col A vector of colors used as a lookup table.
#' @param zero_col Color used for exactly-zero intensities (defaults to transparent).
#' @param alpha Global alpha multiplier applied to all colors when \code{alpha < 1}.
#' @param irange Intensity range used to normalize values before mapping to \code{col}.
#' @param threshold Optional length-2 numeric vector. If \code{diff(threshold) > 0},
#' values within \code{[threshold[1], threshold[2]]} are set to transparent.
#'
#' @return If \code{alpha == 1}, returns a character vector/array of colors.
#' If \code{alpha < 1}, returns an array with an added RGBA channel (last dimension length 4).
#'
#' @importFrom grDevices col2rgb gray heat.colors
#' @export
mapToColors <- function (imslice, col = heat.colors(128, alpha = 1), zero_col = "#00000000",
                         alpha = 1, irange = range(imslice), threshold = c(0, 0)) {

  if (diff(irange) < 0) {
    cli::cli_abort("{.arg irange} must be non-decreasing (diff >= 0), got diff = {diff(irange)}.")
  }
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
