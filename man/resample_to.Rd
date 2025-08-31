\name{resample_to}
\alias{resample_to}
\title{Resample an image with readable method names}
\description{
Convenience wrapper around \code{\link{resample}} that maps
human-friendly method names to the integer \code{interpolation}
argument and (eventually) switches engines.
}
\usage{
resample_to(source, target, method = c("nearest","linear","cubic"),
            engine = c("internal"), ...)
}
\arguments{
\item{source}{A \code{NeuroVol}}
\item{target}{A \code{NeuroVol} or \code{NeuroSpace}}
\item{method}{Interpolation method: \code{"nearest"}, \code{"linear"}, or \code{"cubic"}.}
\item{engine}{Resampling engine. Currently only \code{"internal"} is supported.}
\item{...}{Reserved for future options}
}
\value{A \code{NeuroVol} in the target space.}
\seealso{\code{\link{resample}}}
\examples{
img <- read_vol(system.file("extdata","global_mask_v4.nii", package="neuroim2"))
sp  <- space(img); sp2 <- NeuroSpace(sp@dim*2, sp@spacing/2, origin=sp@origin, trans=trans(img))
r1  <- resample_to(img, sp2, method = "linear")
}

