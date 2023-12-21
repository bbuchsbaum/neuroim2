
#' @keywords internal
#' @noRd
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

}

#' Resample a NeuroVol object
#'
#' This method resamples a NeuroVol object (\code{source}) to match the dimensions and orientation of a target NeuroVol object (\code{target}).
#'
#' @param source A NeuroVol object representing the source volume to be resampled.
#' @param target A NeuroVol object representing the target volume to match the dimensions and orientation of the source volume.
#' @param interpolation A single integer specifying the type of interpolation to be applied to the final resampled image. May be 0 (nearest neighbor), 1 (trilinear), or 3 (cubic spline). No other values are valid.
#' @import RNifti
#' @importFrom RNiftyReg buildAffine applyTransform
#' @return A NeuroVol object with the resampled source volume matched to the target volume's dimensions and orientation.
#' @export
setMethod(f="resample", signature=signature("NeuroVol", "NeuroVol"),
          def=function(source, target, interpolation=3L) {
            hd_target <- as_nifti_header(target, file_name="target.nii")
            hd_source <- as_nifti_header(source, file_name="source.nii")

            hdt <- convert_hd(hd_target)
            hds <- convert_hd(hd_source)

            src <- RNifti::asNifti(as.array(source), reference=hds)
            targ <- RNifti::asNifti(as.array(target), reference=hdt)

            #browser()
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
#' @export
setMethod(f="resample", signature=signature("NeuroVol", "NeuroSpace"),
          def=function(source, target, interpolation=3L) {
            targ <- NeuroVol(array(0, dim(target)), space=target)
            callGeneric(source, targ, interpolation)
          })
