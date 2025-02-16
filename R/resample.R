#' @include all_class.R
NULL

#' Resampling Methods for Neuroimaging Objects
#'
#' @name neuro-resample
#' @description Methods for resampling neuroimaging objects to different spaces and dimensions
NULL

#' @keywords internal
#' @noRd
#' @importFrom RNifti niftiHeader
convert_hd <- function(hd) {
  hdr_targ <- RNifti::niftiHeader(list(
    pixdim=hd$pixdim,
    dim=hd$dimensions,
    qform_code=hd$qform_code,
    quatern_b=hd$quaternion[1],
    quatern_c=hd$quaternion[2],
    quatern_d=hd$quaternion[3],
    qoffset_x=hd$qoffset[1],
    qoffset_y=hd$qoffset[2],
    qoffset_z=hd$qoffset[3]))
  hdr_targ
}

#' Resample a neuroimaging volume
#'
#' @rdname resample-methods
#' @param source The source volume to resample
#' @param target The target space to resample to
#' @return A resampled volume in the target space
#'
#' @examples
#' \dontrun{
#' # Create source and target volumes
#' src_vol <- read_vol(system.file("extdata", "example1.nii", package="neuroim2"))
#' targ_vol <- read_vol(system.file("extdata", "example2.nii", package="neuroim2"))
#'
#' # Resample source to match target
#' resampled <- resample(src_vol, targ_vol, interpolation=1)
#' }
#'
#' @seealso \code{\link{NeuroVol}} for the base volume class
#' @importFrom RNifti asNifti
#' @importFrom RNiftyReg buildAffine applyTransform
#' @export
setMethod(f="resample", signature=signature("NeuroVol", "NeuroVol"),
          def=function(source, target, interpolation=3L) {
            hd_target <- as_nifti_header(target, file_name="target.nii")
            hd_source <- as_nifti_header(source, file_name="source.nii")

            hdt <- convert_hd(hd_target)
            hds <- convert_hd(hd_source)

            src <- RNifti::asNifti(as.array(source), reference=hds)
            targ <- RNifti::asNifti(as.array(target), reference=hdt)

            trans <- RNiftyReg::buildAffine(source=src, target=targ)
            out <- RNiftyReg::applyTransform(trans,src, interpolation=interpolation)

            NeuroVol(unclass(out), space(target))
          })

#' Resample a NeuroVol object to match a NeuroSpace object
#'
#' This method resamples a NeuroVol object (\code{source}) to match the dimensions and orientation of a NeuroSpace object (\code{target}).
#'
#' @param source A NeuroVol object representing the source volume to be resampled.
#' @param target A NeuroSpace object representing the target space to match the dimensions and orientation of the source volume.
#' @param interpolation A single integer specifying the type of interpolation to be applied to the final resampled image. May be 0 (nearest neighbor), 1 (trilinear), or 3 (cubic spline). No other values are valid.
#' @import RNifti
#' @importFrom RNiftyReg buildAffine applyTransform
#' @return A NeuroVol object with the resampled source volume matched to the target space's dimensions and orientation.
#' @rdname resample-methods
#' @export
setMethod(f="resample", signature=signature("NeuroVol", "NeuroSpace"),
          def=function(source, target, interpolation=3L) {
            targ <- NeuroVol(array(0, dim(target)), space=target)
            callGeneric(source, targ, interpolation)
          })
