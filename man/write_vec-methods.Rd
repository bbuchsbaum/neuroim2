% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/latent_neurovec.R,
%   R/neurovec.R
\name{write_vec}
\alias{write_vec}
\alias{write_vec,LatentNeuroVec,character,missing,missing-method}
\alias{write_vec,ROIVec,character,missing,missing-method}
\alias{write_vec,NeuroVec,character,missing,missing-method}
\alias{write_vec,NeuroVec,character,character,missing-method}
\alias{write_vec,NeuroVec,character,missing,character-method}
\alias{write_vec,NeuroVec,character,missing,character,ANY-method}
\title{Write a 4d image vector to disk}
\usage{
write_vec(x, file_name, format, data_type, ...)

\S4method{write_vec}{LatentNeuroVec,character,missing,missing}(x, file_name, nbit = FALSE, compression = 9, chunk_dim = NULL)

\S4method{write_vec}{ROIVec,character,missing,missing}(x, file_name)

\S4method{write_vec}{NeuroVec,character,missing,missing}(x, file_name)

\S4method{write_vec}{NeuroVec,character,character,missing}(
  x,
  file_name,
  format,
  nbit = FALSE,
  compression = 5,
  chunk_dim = c(10, 10, 10, dim(x)[4])
)

\S4method{write_vec}{NeuroVec,character,missing,character}(x, file_name, data_type)
}
\arguments{
\item{x}{an image object, typically a \code{NeuroVec} instance.}

\item{file_name}{output file name.}

\item{format}{file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.}

\item{data_type}{the numeric data type. If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE".
Otherwise output format will be inferred from R the datatype of the image.}

\item{...}{extra args}

\item{nbit}{set nbit compression}

\item{compression}{compression level 1 to 9}

\item{chunk_dim}{the dimensions of each chunk}
}
\description{
Write a 4d image vector to disk
}
\examples{

bvec <- NeuroVec(array(0, c(10,10,10,10)), NeuroSpace(c(10,10,10,10), c(1,1,1)))
\dontrun{
write_vec(bvol, "out.nii")
write_vec(bvol, "out.nii.gz")
write_vec(bvec, "out.nii")
write_vec(bvec, "out.nii.gz")
}
}
