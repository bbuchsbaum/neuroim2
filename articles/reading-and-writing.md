# Reading and writing images

``` r

library(neuroim2)
```

Every analysis starts and ends at a file. This article covers the
readers and writers, how to inspect a header before committing to
loading anything, and the one header problem that produces wrong results
without producing an error.

``` r

mask_file <- system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")
series_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
anat_file <- system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2")
```

## Choosing a reader

| Function | Returns | Use for |
|:---|:---|:---|
| [`read_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md) | `NeuroVol` | a 3D image, or one volume of a 4D one |
| [`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md) | `NeuroVec` | a 4D series |
| [`read_image()`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md) | whichever fits | unknown dimensionality |
| [`read_vol_list()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol_list.md) | `NeuroVec` | many 3D files as one series |
| [`read_header()`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md) | `FileMetaInfo` | metadata without the data |

``` r

vol <- read_vol(mask_file)
vec <- read_vec(series_file)

dim(vol)
#> [1] 64 64 25
dim(vec)
#> [1] 64 64 25  4
```

Compression is handled by extension: `.nii` and `.nii.gz` both work on
read and on write.
[`read_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
on a 4D file returns its first volume, and `index` picks a different
one:

``` r

dim(read_vol(series_file, index = 3))
#> [1] 64 64 25
```

Use
[`read_image()`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md)
when you do not know what a file holds:

``` r

class(read_image(series_file))[1]
#> [1] "DenseNeuroVec"
```

Reading several files gives you a 4D object either way, but not the same
one.
[`read_vol_list()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol_list.md)
builds a single contiguous series;
[`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
on a vector of paths keeps the runs as separable segments:

``` r

class(read_vol_list(c(mask_file, mask_file)))[1]
#> [1] "DenseNeuroVec"
class(read_vec(c(series_file, series_file)))[1]
#> [1] "NeuroVecSeq"

dim(read_vec(c(series_file, series_file)))
#> [1] 64 64 25  8
```

Neither returns a list, despite the name —
[`length()`](https://rdrr.io/r/base/length.html) on either counts
volumes.

## Looking before you load

[`read_header()`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md)
reads metadata only. On a large file this costs nothing and tells you
what you are about to commit to:

``` r

hdr <- read_header(mask_file)

dim(hdr)
#> [1] 64 64 25
hdr@data_type
#> [1] "UBYTE"
hdr@spacing
#> [1] 3.5 3.5 3.7
```

[`dim()`](https://rdrr.io/r/base/dim.html) and
[`trans()`](https://bbuchsbaum.github.io/neuroim2/reference/trans-methods.md)
work on a header, which is the quickest way to see the affine neuroim2
will actually use without loading a byte of data:

``` r

trans(hdr)
#>      [,1] [,2] [,3]    [,4]
#> [1,] -3.5  0.0  0.0  112.00
#> [2,]  0.0  3.5  0.0 -108.50
#> [3,]  0.0  0.0  3.7  -46.25
#> [4,]  0.0  0.0  0.0    1.00
```

[`spacing()`](https://bbuchsbaum.github.io/neuroim2/reference/spacing-methods.md),
[`origin()`](https://bbuchsbaum.github.io/neuroim2/reference/origin-methods.md)
and
[`space()`](https://bbuchsbaum.github.io/neuroim2/reference/space-methods.md)
are not defined for headers — read the slots (`hdr@spacing`,
`hdr@origin`) instead. The full NIfTI header is a plain list in the
`header` slot, so anything the slots omit is still reachable:

``` r

hdr@header$datatype
#> [1] 2
hdr@header$encoding
#> [1] "gzip"
hdr@header$vox_offset
#> [1] 3040
```

## The two affines

A NIfTI file can carry **two** spatial transforms: the `qform` (a rigid
transform stored as a quaternion) and the `sform` (a general affine).
Each has a code saying whether it is set. The standard says use the
sform when `sform_code > 0`, otherwise the qform.

On a well-formed file the two agree:

``` r

c(qform_code = hdr@header$qform_code, sform_code = hdr@header$sform_code)
#> qform_code sform_code 
#>          1          1

all.equal(hdr@header$qform, hdr@header$sform)
#> [1] TRUE
```

They do not always agree. The anatomical image shipped with this package
is a real example — its qform describes 4 mm voxels and its sform
describes 1 mm ones:

``` r

anat_hdr <- read_header(anat_file)

anat_hdr@header$pixdim[2:4]
#> [1] 4.02083 4.01754 4.02083
anat_hdr@header$qform
#>         [,1]    [,2]    [,3] [,4]
#> [1,] 4.02083 0.00000 0.00000  -96
#> [2,] 0.00000 4.01754 0.00000 -132
#> [3,] 0.00000 0.00000 4.02083  -78
#> [4,] 0.00000 0.00000 0.00000    1
anat_hdr@header$sform
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0  -96
#> [2,]    0    1    0 -132
#> [3,]    0    0    1  -78
#> [4,]    0    0    0    1
```

Both codes are set, so the sform wins and the image reports 1 mm voxels
even though `pixdim` says otherwise:

``` r

spacing(read_vol(anat_file))
#> [1] 1 1 1
```

Nothing errors. The image loads and plots correctly, and every
**distance, volume and ROI radius** derived from it is wrong by a factor
of about four. The coordinates themselves are not scaled by four — both
affines share a translation, so they agree exactly at voxel `(1, 1, 1)`
and diverge by a growing offset that reaches roughly 140 mm at the far
corner:

``` r

far <- matrix(dim(read_vol(anat_file)), nrow = 1)

sform_xyz <- grid_to_coord(space(read_vol(anat_file)), far)
qform_xyz <- (anat_hdr@header$qform %*% c(far - 1, 1))[1:3]

rbind(sform = as.vector(sform_xyz), qform = qform_xyz)
#>            [,1]      [,2]    [,3]
#> sform -49.00000 -76.00000 -31.000
#> qform  92.97902  92.98224 110.979
```

**This is the failure mode to know about**, because no amount of care
downstream will catch it. Check when a file first arrives:

``` r

affines_agree <- function(file) {
  h <- read_header(file)
  qc <- h@header$qform_code
  sc <- h@header$sform_code
  if (is.null(qc) || is.null(sc) || qc <= 0 || sc <= 0) return(NA)
  isTRUE(all.equal(h@header$qform, h@header$sform, tolerance = 1e-4))
}

c(mask = affines_agree(mask_file), anat = affines_agree(anat_file))
#>  mask  anat 
#>  TRUE FALSE
```

The [`is.null()`](https://rdrr.io/r/base/NULL.html) guard matters: not
every format carries these fields, and returning `NA` for “cannot tell”
is better than erroring on an AFNI file.

If a file is affected and you trust its `pixdim`, rebuild the geometry:

``` r

bad <- read_vol(anat_file)
fixed <- NeuroVol(
  as.array(bad),
  NeuroSpace(dim(bad),
    spacing = anat_hdr@header$pixdim[2:4],
    origin = anat_hdr@header$qform[1:3, 4],
    axes = space(bad)@axes
  )
)

spacing(fixed)
#> [1] 4.02083 4.01754 4.02083
```

## Writing

[`write_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.md)
and
[`write_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vec-methods.md)
take an object and a path. The extension decides the container and the
compression; an unrecognised or absent extension gives uncompressed
NIfTI.

``` r

out_nii <- tempfile(fileext = ".nii")
out_gz <- tempfile(fileext = ".nii.gz")

write_vol(vol, out_nii)
write_vol(vol, out_gz)

c(plain = file.size(out_nii), gzipped = file.size(out_gz))
#>   plain gzipped 
#>  102752    2157
```

That plain file is 400 kB for a binary mask, because
[`write_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.md)
promotes to `FLOAT` by default. Pass `data_type` when you know better:

``` r

out_byte <- tempfile(fileext = ".nii")
write_vol(vol, out_byte, data_type = "UBYTE")

c(float = file.size(out_nii), ubyte = file.size(out_byte))
#>  float  ubyte 
#> 102752 102752
read_header(out_byte)@data_type
#> [1] "UBYTE"
```

Valid types are `UNKNOWN`, `BINARY`, `UBYTE`, `SHORT`, `INT`, `FLOAT`
and `DOUBLE`. A round trip preserves both data and geometry either way:

``` r

back <- read_vol(out_gz)

all.equal(as.array(back), as.array(vol))
#> [1] TRUE
identical(trans(space(back)), trans(space(vol)))
#> [1] TRUE
```

4D objects work the same way:

``` r

out_vec <- tempfile(fileext = ".nii.gz")
write_vec(sub_vector(vec, 1:2), out_vec)

dim(read_vec(out_vec))
#> [1] 64 64 25  2
```

### A precision surprise that is not about files

Affines round to seven significant figures, and spacing and origin to
six — but this happens in
[`NeuroSpace()`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace.md),
when the object is built, not when it is written:

``` r

sp <- NeuroSpace(c(4L, 4L, 4L), spacing = c(1 / 3, 0.123456789012, pi))

spacing(sp)
#> [1] 0.333333 0.123457 3.141590
```

No file was involved. Because both sides of a round trip get the same
rounding, [`identical()`](https://rdrr.io/r/base/identical.html) on a
re-read affine does succeed, as above. What you cannot rely on is an
affine you constructed matching the full-precision numbers you passed
in.

## AFNI files

AFNI `HEAD`/`BRIK` pairs are read by the same functions — pass the
`.HEAD` file and the matching `.BRIK` is found beside it.
[`read_header()`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md)
returns an `AFNIMetaInfo` whose `header` slot holds the AFNI attribute
list:

``` r

vol <- read_vol("subject_anat+orig.HEAD")
names(read_header("subject_anat+orig.HEAD")@header)
```

## Where to go next

- [`vignette("spaces-and-coordinates")`](https://bbuchsbaum.github.io/neuroim2/articles/spaces-and-coordinates.md)
  — what the affine you just read means
- [`vignette("large-data")`](https://bbuchsbaum.github.io/neuroim2/articles/large-data.md)
  — reading files too big to hold in memory
- [`?read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md),
  [`?read_header`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md),
  [`?write_vol`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.md)
