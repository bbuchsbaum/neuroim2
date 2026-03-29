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

  The name(s) of the file(s) to load. If multiple files are specified,
  they are loaded and concatenated along the time dimension.

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

An
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
object representing the loaded volume(s).

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
