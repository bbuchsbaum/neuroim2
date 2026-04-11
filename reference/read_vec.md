# read_vec

Loads a neuroimaging volume from one or more files, with support for
various input formats and memory management strategies.

## Usage

``` r
read_vec(
  file_name,
  indices = NULL,
  mask = NULL,
  mode = c("normal", "mmap", "bigvec", "filebacked")
)
```

## Arguments

- file_name:

  A character vector of one or more file paths to load. A 3D file is
  promoted to a 4D `NeuroVec` with a single time point (see Details).
  When multiple paths are supplied the result is always a
  [`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md)
  (which itself extends
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)),
  regardless of whether the individual files are 3D, 4D, or a mix of
  both.

- indices:

  The indices of the sub-volumes to load (e.g. if the file is
  4-dimensional). Only supported in "normal" mode.

- mask:

  A logical mask defining which spatial elements to load. Required for
  "bigvec" mode and optional for other modes.

- mode:

  The IO mode which is one of: \* "normal": Standard in-memory loading
  \* "mmap": Memory-mapped access (more memory efficient) \* "bigvec":
  Optimized for large datasets with masking \* "filebacked": File-backed
  storage with on-demand loading

## Value

The return type depends on how many files are supplied:

- **Single 3D file** — a
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  with `dim(x)[4] == 1` (concrete class depends on `mode`: e.g.
  [`DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md)
  for `"normal"`,
  [`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
  for `"mmap"`,
  [`BigNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/BigNeuroVec-class.md)
  for `"bigvec"`,
  [`FileBackedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/FileBackedNeuroVec-class.md)
  for `"filebacked"`).

- **Single 4D file** — a
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  of the same concrete class as above, with `dim(x)[4]` equal to the 4th
  dimension of the file (or `length(indices)` when `indices` is
  supplied).

- **Multiple files (any mix of 3D and 4D)** — a
  [`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md)
  wrapping one `NeuroVec` per input file in the order given. Because
  `NeuroVecSeq` extends
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
  it can be used wherever a `NeuroVec` is accepted, but its underlying
  storage is segmented rather than a single contiguous 4D array. For
  example, `read_vec(c(vol, vec, vec, vol))` returns a 4-element
  `NeuroVecSeq` whose segments have time lengths `c(1, T2, T3, 1)`.

## Details

This function supports multiple file formats: \* .nii: Standard NIfTI
format \* .nii.gz: Compressed NIfTI (not supported in mmap mode)

Memory management modes: \* "normal": Loads entire dataset into memory.
Best for smaller datasets or when memory is not a constraint. \* "mmap":
Memory-maps the file, providing efficient access for large files without
loading entirely into memory. Not available for compressed files. \*
"bigvec": Optimized for large datasets where only a subset of voxels are
of interest. Requires a mask to specify which voxels to load. \*
"filebacked": Similar to mmap but with more flexible caching strategies.

**3D inputs:** A path pointing at a 3D image is not rejected. It is
promoted to a 4D `NeuroVec` whose fourth dimension has length 1, so the
return type is always a `NeuroVec`, never a
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md).
If you want a true volume, use
[`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
(or `vec[[1]]`).

**Multiple files:** When `file_name` has length \> 1, each file is
loaded independently and the results are wrapped with
[`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq.md).
No data is copied or re-concatenated into a single dense 4D array; the
returned `NeuroVecSeq` holds the constituent `NeuroVec`s in its `@vecs`
slot and exposes them as a single logical time series. Time-point
lookups (e.g. `x[[i]]`, `sub_vector(x, i)`, `linear_access(x, i)`)
transparently walk across the segments. The per-file time lengths may
differ (e.g. a 3D file contributes 1 time point, a 4D file contributes
`dim(file)[4]`); all spatial dimensions must match.

## Note

\* Memory-mapping (.mmap mode) is not supported for gzipped files \* The
bigvec mode requires a mask to be specified \* When loading multiple
files, they must have compatible dimensions

## Examples

``` r
# Load a single NIfTI file
img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))


# Memory-mapped loading for large files
big_img <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"), mode="mmap")
#> loading /home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii as mmaped file 

# Load masked data for memory efficiency
mask <- as.logical(big_img[[1]])
masked_data <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"),
               mask=mask, mode="bigvec")
#> loading /home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii as bigvec 

```
