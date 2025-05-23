% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spat_filter.R
\name{gaussian_blur}
\alias{gaussian_blur}
\title{Gaussian Blur for Volumetric Images}
\usage{
gaussian_blur(vol, mask, sigma = 2, window = 1)
}
\arguments{
\item{vol}{A \code{\linkS4class{NeuroVol}} object representing the image volume to be smoothed.}

\item{mask}{An optional \code{\linkS4class{LogicalNeuroVol}} object representing the image mask.
This mask defines the region where the blurring is applied. If not provided, the entire volume is processed.}

\item{sigma}{A numeric value specifying the standard deviation of the Gaussian kernel. Default is 2.}

\item{window}{An integer specifying the kernel size. It represents the number of voxels to include
on each side of the center voxel. For example, window=1 results in a 3x3x3 kernel. Default is 1.}
}
\value{
A \code{\linkS4class{NeuroVol}} object representing the smoothed image.
}
\description{
This function applies an isotropic discrete Gaussian kernel to smooth a volumetric image (3D brain MRI data).
The blurring is performed within a specified image mask, with customizable kernel parameters.
}
\details{
The function uses a C++ implementation for efficient Gaussian blurring. The blurring is applied
only to voxels within the specified mask (or the entire volume if no mask is provided).
The kernel size is determined by the 'window' parameter, and its shape by the 'sigma' parameter.
}
\examples{
# Load a sample brain mask
brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Apply Gaussian blurring to the brain volume
blurred_vol <- gaussian_blur(brain_mask, brain_mask, sigma = 2, window = 1)

# View a slice of the original and blurred volumes
image(brain_mask[,,12])
image(blurred_vol[,,12])

}
\references{
Gaussian blur: https://en.wikipedia.org/wiki/Gaussian_blur
}
\seealso{
\code{\link{NeuroVol-class}}, \code{\link{LogicalNeuroVol-class}}, \code{\link{bilateral_filter}}
}
