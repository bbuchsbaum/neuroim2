% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurospace.R
\docType{methods}
\name{origin}
\alias{origin}
\alias{origin,NeuroSpace-method}
\title{Generic getter to extract image origin}
\usage{
origin(x)

\S4method{origin}{NeuroSpace}(x)
}
\arguments{
\item{x}{an object with an origin}
}
\description{
Generic getter to extract image origin
}
\examples{
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
origin(bspace)

}