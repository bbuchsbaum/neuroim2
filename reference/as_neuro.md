# Convert between neuroim2 objects and RNifti / oro.nifti images

`as_neuro()` turns an image held by another R neuroimaging package into
a neuroim2 object, and `as_nifti()` goes the other way. Together they
let neuroim2 take part in a pipeline built on RNifti, RNiftyReg,
oro.nifti, neurobase or fslr without a round trip through a file.

## Usage

``` r
as_neuro(x, ..., label = "")

as_nifti(x, xform_code = 1L)

# S3 method for class 'NeuroObj'
asNifti(x, reference = NULL, ..., xform_code = 1L)
```

## Arguments

- x:

  For `as_neuro()`, a `niftiImage`, `internalImage`, oro.nifti `nifti`
  object, file name or neuroim2 object. For `as_nifti()`, a
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  or
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md).

- ...:

  For `as_neuro()`, further arguments passed to
  [`read_image`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md)
  when `x` is a file name; otherwise unused.

- label:

  For `as_neuro()`, an optional label for the result.

- xform_code:

  The code written to `sform_code` and `qform_code`: 1 (scanner
  anatomical, the default), 2 (aligned to another image), 3 (Talairach)
  or 4 (MNI152). neuroim2 does not record which of these a space is, so
  say so here when it matters downstream.

- reference:

  Ignored; present because
  [`RNifti::asNifti()`](https://rdrr.io/pkg/RNifti/man/asNifti.html) has
  it.

## Value

`as_neuro()` returns a
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
or
[`DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md).
`as_nifti()` returns a `niftiImage`.

## Details

`as_neuro()` accepts

- an RNifti `niftiImage` or `internalImage`, e.g. the result of
  [`RNifti::readNifti()`](https://rdrr.io/pkg/RNifti/man/readNifti.html)
  or of an RNiftyReg registration;

- an oro.nifti `nifti` object, e.g. the result of
  [`oro.nifti::readNIfTI()`](https://rdrr.io/pkg/oro.nifti/man/read_nifti.html),
  `neurobase::readnii()` or an fslr command;

- a file name, which is passed to
  [`read_image`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md);

- a neuroim2 object, which is returned unchanged.

A 3-D image becomes a
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
and a 4-D image a
[`DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md).
Trailing dimensions of extent 1 are dropped first, so a 4-D image with
one volume becomes a volume.

The voxel-to-world affine is taken from
[`RNifti::xform()`](https://rdrr.io/pkg/RNifti/man/xform.html), which
applies the NIfTI rule: the sform if it is set, else the qform, else a
scaling by the voxel sizes. RNifti and neuroim2 share one convention for
it (zero-based voxel indices to RAS+ millimetres), so
`trans(as_neuro(img))` equals `RNifti::xform(img)` to the seven
significant figures a
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
keeps. Data scaling (`scl_slope`, `scl_inter`) has already been applied
by the package that read the image, so values carry over as they are.

`as_nifti()` returns an RNifti `niftiImage` holding the values of `x` as
a dense array. The affine of `x` is written to both the sform and the
qform, and the voxel sizes to `pixdim`. The qform can only describe a
rotation with scaling, so for an affine with shear the qform is the
closest such transform and the sform is the exact one; readers that
follow the NIfTI rule use the sform. Any object with an
[`as.array()`](https://bbuchsbaum.github.io/neuroim2/reference/as.array.md)
method and a
[`space`](https://bbuchsbaum.github.io/neuroim2/reference/space-methods.md)
can be converted, so sparse, clustered, memory-mapped and file-backed
objects are read into memory by the call. A logical volume is stored as
0/1 integers.

Because
[`RNifti::asNifti()`](https://rdrr.io/pkg/RNifti/man/asNifti.html) is an
S3 generic, `asNifti(x)` works on neuroim2 objects as well and gives the
same result as `as_nifti(x)`. To get an oro.nifti object, pass the
result on: `oro.nifti::nii2oro(as_nifti(x))`.

Volume labels, NIfTI extensions and header fields other than the
geometry are not carried across in either direction.

## See also

[`read_image`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md),
[`write_vol`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.md),
[`as_nifti_header`](https://bbuchsbaum.github.io/neuroim2/reference/as_nifti_header.md)

## Examples

``` r
fname <- system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2")

# An image read by RNifti becomes a neuroim2 volume ...
img <- RNifti::readNifti(fname)
vol <- as_neuro(img)
all.equal(trans(vol), RNifti::xform(img), check.attributes = FALSE,
          tolerance = 1e-6)
#> [1] TRUE

# ... and a neuroim2 volume goes back, with its geometry
smoothed <- gaussian_blur(vol, as.mask(vol > 0), sigma = 2)
img2 <- as_nifti(smoothed)
RNifti::pixdim(img2)
#> [1] 4.02083 4.01754 4.02083
```
