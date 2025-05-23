% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neurovol.R
\name{read_vol}
\alias{read_vol}
\title{Load an image volume from a file}
\usage{
read_vol(file_name, index = 1)
}
\arguments{
\item{file_name}{the name of the file to load}

\item{index}{the index of the volume (e.g. if the file is 4-dimensional)}
}
\value{
an instance of the class \code{\linkS4class{DenseNeuroVol}}
}
\description{
Load an image volume from a file
}
\examples{
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
x <- read_vol(fname)
print(dim(x))
space(x)

}
