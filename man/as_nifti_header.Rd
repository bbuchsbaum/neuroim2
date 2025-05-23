% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nifti_io.R
\name{as_nifti_header}
\alias{as_nifti_header}
\title{Construct a Minimal NIfTI-1 Header from a NeuroVol}
\usage{
as_nifti_header(vol, file_name, oneFile = TRUE, data_type = "FLOAT")
}
\arguments{
\item{vol}{A \code{\link[neuroim2]{NeuroVol}} (or 3D array-like) specifying
dimensions, spacing, and affine transform.}

\item{file_name}{A character string for the file name (used within the header
but not necessarily to write data).}

\item{oneFile}{Logical; if \code{TRUE}, sets the NIfTI magic to \code{"n+1"},
implying a single-file format (\code{.nii}). If \code{FALSE}, uses
\code{"ni1"} (header+image).}

\item{data_type}{Character specifying the data representation, e.g. \code{"FLOAT"},
\code{"DOUBLE"}. The internal code picks an integer NIfTI code.}
}
\value{
A \code{list} representing the NIfTI-1 header fields, containing
  elements like \code{dimensions}, \code{pixdim}, \code{datatype}, \code{qform},
  \code{quaternion}, \code{qfac}, etc. This can be passed to other
  functions that write or manipulate the header.
}
\description{
Given a \code{\link[neuroim2]{NeuroVol}} object (or similar), this function
builds a basic NIfTI-1 header structure, populating essential fields such as
\code{dim}, \code{pixdim}, \code{datatype}, the affine transform, and the
quaternion parameters.
}
\details{
This is a convenience function that calls \code{\link{createNIfTIHeader}}
first, then updates the fields (dimensions, \code{pixdim}, orientation, etc.)
based on the \code{vol} argument. The voxel offset is set to 352 bytes, and
the quaternion is derived from the transform matrix via
\code{\link{matrixToQuatern}}.

Note: This function primarily sets up a minimal header suitable for writing
standard single-file NIfTI-1. If you need a more comprehensive or advanced
usage, consider manually editing the returned list.
}
\seealso{
\code{\link{createNIfTIHeader}} for the base constructor of an empty NIfTI header.
}
