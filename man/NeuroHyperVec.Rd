% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neurohypervec.R
\name{NeuroHyperVec}
\alias{NeuroHyperVec}
\title{Create a NeuroHyperVec instance}
\usage{
NeuroHyperVec(...)
}
\arguments{
\item{...}{One or more instances of type NeuroVec}
}
\value{
A NeuroHyperVec object
}
\description{
Create a NeuroHyperVec instance
}
\examples{
space <- NeuroSpace(dim = c(64, 64, 32), origin = c(0, 0, 0), spacing = c(3, 3, 4))
vec1 <- NeuroVec(data = array(rnorm(64*64*32*100), dim = c(64, 64, 32, 100)), space = space)
vec2 <- NeuroVec(data = array(rnorm(64*64*32*100), dim = c(64, 64, 32, 100)), space = space)
hyper_vec <- NeuroHyperVec(vec1, vec2)
}
