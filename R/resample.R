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

#' @keywords internal
#' @noRd
.resample_clustered_neurovol <- function(source, target_space, interpolation=0L) {
  if (!identical(interpolation, 0L)) {
    warning("ClusteredNeuroVol resampling uses nearest-neighbor interpolation; forcing interpolation = 0.",
            call. = FALSE)
  }

  # resample labels using nearest neighbor to keep integer ids intact
  label_vol <- NeuroVol(as.numeric(as.array(as.dense(source))), space(source))
  resampled_labels <- resample(label_vol, target_space, interpolation=0L)
  label_array <- as.array(resampled_labels)

  # resample mask to ensure we stay within the original support
  mask_vol <- NeuroVol(as.numeric(as.array(mask(source))), space(source))
  resampled_mask <- resample(mask_vol, target_space, interpolation=0L)
  mask_array <- as.array(resampled_mask)

  label_array[is.na(label_array)] <- 0
  mask_array[is.na(mask_array)] <- 0

  keep_mask <- (mask_array != 0) & (label_array != 0)
  target_mask <- LogicalNeuroVol(keep_mask, target_space)

  if (!any(keep_mask)) {
    return(ClusteredNeuroVol(target_mask, integer(0), label_map=list()))
  }

  clusters <- as.integer(round(label_array[keep_mask]))
  ids_present <- sort(unique(clusters))

  label_map <- source@label_map
  if (length(label_map) > 0) {
    keep <- vapply(label_map, function(v) any(unlist(v) %in% ids_present),
                   logical(1))
    label_map <- label_map[keep]
  } else {
    label_map <- NULL
  }

  ClusteredNeuroVol(target_mask, clusters, label_map=label_map)
}

#' Resample a neuroimaging volume
#'
#' @rdname resample-methods
#' @param source The source volume to resample
#' @param target The target space to resample to
#'
#' @examples
#' # Create source and target volumes
#' src_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#' targ_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
#'
#' # Resample source to match target
#' resampled <- resample(src_vol, targ_vol, interpolation=1)
#' 
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
#' @rdname resample-methods
#' @export
setMethod(f="resample", signature=signature("NeuroVol", "NeuroSpace"),
          def=function(source, target, interpolation=3L) {
            targ <- NeuroVol(array(0, dim(target)), space=target)
            callGeneric(source, targ, interpolation)
          })

#' Resample a ClusteredNeuroVol using nearest-neighbor interpolation
#'
#' This method preserves discrete cluster labels and label mappings when
#' resampling clustered volumes to a new space.
#'
#' @rdname resample-methods
#' @export
setMethod(f="resample", signature=signature("ClusteredNeuroVol", "NeuroSpace"),
          def=function(source, target, interpolation=0L) {
            .resample_clustered_neurovol(source, target, interpolation)
          })

#' @rdname resample-methods
#' @export
setMethod(f="resample", signature=signature("ClusteredNeuroVol", "NeuroVol"),
          def=function(source, target, interpolation=0L) {
            .resample_clustered_neurovol(source, space(target), interpolation)
          })
