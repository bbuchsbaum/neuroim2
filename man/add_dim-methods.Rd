% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurospace.R
\name{add_dim}
\alias{add_dim}
\alias{add_dim,NeuroSpace,numeric-method}
\title{Add a Dimension to an Object}
\usage{
add_dim(x, n)

\S4method{add_dim}{NeuroSpace,numeric}(x, n)
}
\arguments{
\item{x}{A dimensioned object, such as a matrix, an array, or a NeuroSpace object.}

\item{n}{An integer representing the size of the dimension to add.}
}
\value{
An updated dimensioned object with the new dimension added.
}
\description{
This function adds a new dimension to a given object, such as a matrix or an array.
}
\examples{
# Create a NeuroSpace object
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Add a new dimension with size 10
x1 <- add_dim(x, 10)

# Check the new dimension
ndim(x1) == 4
dim(x1)[4] == 10

}
