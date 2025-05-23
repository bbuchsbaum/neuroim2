% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neurovecseq.R
\name{[[,NeuroVecSeq,numeric-method}
\alias{[[,NeuroVecSeq,numeric-method}
\title{Extract Element from NeuroVecSeq}
\usage{
\S4method{[[}{NeuroVecSeq,numeric}(x, i)
}
\arguments{
\item{x}{A NeuroVecSeq object}

\item{i}{Numeric index specifying the time point to extract}
}
\value{
A NeuroVol object representing the extracted volume
}
\description{
Extracts a single volume from a NeuroVecSeq object at the specified time point.
}
