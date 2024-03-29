% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/common.R
\name{read_columns}
\alias{read_columns}
\alias{read_columns,ColumnReader,integer-method}
\title{Read a set of column vector from an input source (e.g. \code{ColumnReader})}
\usage{
read_columns(x, column_indices)

\S4method{read_columns}{ColumnReader,integer}(x, column_indices)
}
\arguments{
\item{x}{the input channel}

\item{column_indices}{the column indices}
}
\value{
a \code{matrix} consisting of the requested column vectors
}
\description{
Read a set of column vector from an input source (e.g. \code{ColumnReader})
}
