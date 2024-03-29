% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R
\docType{class}
\name{FileBackedNeuroVec-class}
\alias{FileBackedNeuroVec-class}
\title{FileBackedNeuroVec}
\description{
A class representing a four-dimensional brain image that is read into memory "on demand" using memory-mapped file access. FileBackedNeuroVec objects enable efficient handling of large-scale brain imaging data by loading only the required portions of the data into memory when needed.
}
\section{Slots}{

\describe{
\item{\code{meta}}{An instance of class \code{\linkS4class{FileMetaInfo}} containing the file meta information, such as the file path, format, and other associated metadata.}
}}

