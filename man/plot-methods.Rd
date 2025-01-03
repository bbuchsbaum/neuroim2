% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neuroslice.R, R/neurovol.R
\name{plot,NeuroSlice-method}
\alias{plot,NeuroSlice-method}
\alias{plot,NeuroVol-method}
\title{Visualize 2D Neuroimaging Slice}
\usage{
\S4method{plot}{NeuroSlice}(
  x,
  cmap = gray(seq(0, 1, length.out = 255)),
  irange = range(x, na.rm = TRUE)
)

\S4method{plot}{NeuroVol}(
  x,
  cmap = gray(seq(0, 1, length.out = 255)),
  zlevels = unique(round(seq(1, dim(x)[3], length.out = 6))),
  irange = range(x, na.rm = TRUE),
  thresh = c(0, 0),
  alpha = 1,
  bgvol = NULL,
  bgcmap = gray(seq(0, 1, length.out = 255))
)
}
\arguments{
\item{x}{the object to display}

\item{cmap}{a color map consisting of a vector of colors in hex format (e.g. \code{gray(n=255)})}

\item{irange}{the intensity range indicating the low and high values of the color scale.}

\item{zlevels}{the series of slice indices to display.}

\item{thresh}{a 2-element vector indicating the lower and upper transparency thresholds.}

\item{alpha}{the level of alpha transparency}

\item{bgvol}{a background volume that serves as an image underlay (currently ignored).}

\item{bgcmap}{a color map for backround layer consisting of a vector of colors in hex format (e.g. \code{gray(n=255)})}

\item{...}{Additional arguments passed to plotting methods}
}
\value{
A \code{ggplot2} object
}
\description{
Creates a 2D visualization of a \code{NeuroSlice} object using \code{ggplot2}.
}
\details{
Plot a NeuroSlice Object


The plot method uses \code{ggplot2} to create a raster visualization of the slice data.
The intensity values are mapped to colors using the specified colormap and range.
}
\examples{
# Create example slice
slice_space <- NeuroSpace(c(100, 100))
slice_data <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(slice_data, slice_space)

\dontrun{
# Basic plot
plot(slice)

# Custom colormap
library(viridis)
plot(slice, cmap = viridis(255))

# Custom intensity range
plot(slice, irange = c(-2, 2))
}


dat <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(dat, NeuroSpace(c(100,100)))
#plot(slice)
}
