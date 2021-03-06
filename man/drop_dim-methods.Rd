% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/axis.R, R/neurospace.R
\name{drop_dim}
\alias{drop_dim}
\alias{drop_dim,AxisSet2D,numeric-method}
\alias{drop_dim,AxisSet2D,missing-method}
\alias{drop_dim,AxisSet3D,numeric-method}
\alias{drop_dim,AxisSet3D,missing-method}
\alias{drop_dim,NeuroSpace,numeric-method}
\alias{drop_dim,NeuroSpace,missing-method}
\title{drop a dimension from an object}
\usage{
drop_dim(x, dimnum)

\S4method{drop_dim}{AxisSet2D,numeric}(x, dimnum)

\S4method{drop_dim}{AxisSet2D,missing}(x, dimnum)

\S4method{drop_dim}{AxisSet3D,numeric}(x, dimnum)

\S4method{drop_dim}{AxisSet3D,missing}(x, dimnum)

\S4method{drop_dim}{NeuroSpace,numeric}(x, dimnum)

\S4method{drop_dim}{NeuroSpace,missing}(x)
}
\arguments{
\item{x}{a dimensioned object}

\item{dimnum}{the index of the dimension to drop}
}
\description{
drop a dimension from an object
}
\examples{
x = NeuroSpace(c(10,10,10), c(1,1,1))
x1 <- drop_dim(x)
ndim(x1) == 2
dim(x1)[2] == 10
}
