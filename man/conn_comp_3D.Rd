% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conncomp.R
\name{conn_comp_3D}
\alias{conn_comp_3D}
\title{Extract connected components from a 3D mask}
\usage{
conn_comp_3D(mask, connect = c("26-connect", "18-connect", "6-connect"))
}
\arguments{
\item{mask}{a 3D binary array}

\item{connect}{the connectiivty constraint: "6-connect", "18-connect", or "26-connect"}
}
\value{
a two-element list of the connected components (cluster \code{index} and cluster \code{size})
The first element \code{index} is a 3D array containing the cluster index of the connected component for each voxel.
The second element \code{size} is a 3D array consisting of the size of the connected component inhabited by each voxel.
}
\description{
Extract connected components from a 3D mask
}
\examples{

dat <- array(as.logical(rnorm(10*10*10)>.5), c(10, 10, 10))
res1 <- conn_comp_3D(dat, connect="6-connect")
res2 <- conn_comp_3D(dat, connect="18-connect")
res3 <- conn_comp_3D(dat, connect="26-connect")

}
