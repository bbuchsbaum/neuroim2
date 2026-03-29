# Access NIfTI Header Information

Retrieves header metadata from neuroimaging objects or files. Returns a
structured list with commonly needed fields like sform/qform matrices,
TR, intent codes, and data scaling parameters.

For low-level access to all raw NIfTI header fields, use the `$raw`
element of the returned list.

## Usage

``` r
header(x)

# S4 method for class 'FileMetaInfo'
header(x)

# S4 method for class 'character'
header(x)
```

## Arguments

- x:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md),
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
  [`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md),
  or a character file path.

## Value

A list of class `"NeuroHeader"` with elements:

- dim:

  Integer dimensions.

- pixdim:

  Voxel sizes including TR.

- spacing:

  Spatial voxel sizes (first 3 pixdim values).

- origin:

  Coordinate origin.

- trans:

  4x4 affine transform.

- qform:

  List with `matrix` (4x4) and `code` (integer).

- sform:

  List with `matrix` (4x4) and `code` (integer).

- intent_code:

  NIfTI intent code.

- intent_name:

  NIfTI intent name string.

- descrip:

  Description string.

- data_type:

  Storage data type label.

- bitpix:

  Bits per pixel.

- scl_slope:

  Data scaling slope.

- scl_inter:

  Data scaling intercept.

- cal_min:

  Display intensity minimum.

- cal_max:

  Display intensity maximum.

- TR:

  Repetition time (4th pixdim), or `NA` if not 4D.

- raw:

  The complete raw header list (all NIfTI fields).

## Examples

``` r
# \donttest{
f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
h <- header(f)
h$dim
#> [1] 64 64 25  4
h$TR
#> [1] 0
h$sform
#> $matrix
#>      [,1] [,2] [,3]   [,4]
#> [1,] -3.5  0.0  0.0  112.0
#> [2,]  0.0  3.5  0.0 -108.0
#> [3,]  0.0  0.0  3.7  -46.2
#> [4,]  0.0  0.0  0.0    1.0
#> 
#> $code
#> [1] 1
#> 
h$descrip
#> [1] NA
# }
```
