# Repair the sform of inst/extdata/mni_downsampled.nii.gz (#50).
#
# The file is the MNI152 1 mm template (193 x 229 x 193) downsampled to
# 48 x 57 x 48, i.e. ~4.02 mm voxels. Its pixdim and qform say so, but its
# sform (sform_code = 1, which takes precedence) was the identity plus the
# translation, so readers reported 1 mm voxels. This copies the qform into the
# sform. Voxel data and datatype are unchanged.
#
# Run from the package root: Rscript tools/fix-mni-downsampled-sform.R

path <- file.path("inst", "extdata", "mni_downsampled.nii.gz")
img <- RNifti::readNifti(path, internal = TRUE)

q <- RNifti::xform(img, useQuaternionFirst = TRUE)
stopifnot(attr(q, "code") > 0)
RNifti::sform(img) <- structure(q[], code = 1L)

tmp <- tempfile(fileext = ".nii.gz")
RNifti::writeNifti(img, tmp, datatype = "float")

before <- RNifti::readNifti(path)
after <- RNifti::readNifti(tmp)
stopifnot(
  identical(dim(before), dim(after)),
  isTRUE(all.equal(as.vector(before), as.vector(after))),
  isTRUE(all.equal(
    RNifti::xform(after, useQuaternionFirst = FALSE)[],
    RNifti::xform(after, useQuaternionFirst = TRUE)[],
    check.attributes = FALSE
  ))
)
file.copy(tmp, path, overwrite = TRUE)
