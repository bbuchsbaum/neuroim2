convert_hd <- function(hd) {
  hdr_targ <- RNifti::niftiHeader(list(
    pixdim=hd$pixdim,
    dimensions=hd$dimensions,
    qform_code=hd$qform_code,
    quatern_b=hd$quaternion[1],
    quatern_c=hd$quaternion[2],
    quatern_d=hd$quaternion[3],
    qoffset_x=hd$qoffset[1],
    qoffset_y=hd$qoffset[2],
    qoffset_z=hd$qoffset[3]))

}

#' @export
#' @import RNifti
#' @importFrom RNiftyReg buildAffine applyTransform
#' @rdname centroids-methods
#' @param interpolation a single integer specifying the type of interpolation to be applied to the
#' final resampled image. May be 0 (nearest neighbour), 1 (trilinear) or 3 (cubic spline).
#' No other values are valid.
setMethod(f="resample", signature=signature("NeuroVol", "NeuroVol"),
          def=function(source, target, interpolation=3L) {
            hd_target <- as_nifti_header(target, file_name="target.nii")
            hd_source <- as_nifti_header(source, file_name="source.nii")

            hdt <- convert_hd(hd_target)
            hds <- convert_hd(hd_source)

            src <- RNifti::asNifti(as.array(source), reference=hds)
            targ <- RNifti::asNifti(as.array(target), reference=hdt)

            trans <- RNiftyReg::buildAffine(source=src, target=targ)
            out <- RNiftyReg::applyTransform(trans,src, interpolation=interpolation, nearest=nearest)

            NeuroVol(unclass(out), space(target))

          })
