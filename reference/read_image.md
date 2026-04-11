# read_image

Convenience wrapper that inspects the file header(s) and dispatches to
the appropriate specialized reader:
[`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
for a single 3D file,
[`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
for 4D data or for any multi-file input, or
[`read_hyper_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_hyper_vec.md)
for 5D data.

## Usage

``` r
read_image(
  file_name,
  type = c("auto", "vol", "vec", "hyper"),
  index = 1,
  indices = NULL,
  mask = NULL,
  mode = c("normal", "mmap", "bigvec", "filebacked")
)
```

## Arguments

- file_name:

  Character vector of one or more file paths.

- type:

  One of `"auto"`, `"vol"`, `"vec"`, or `"hyper"` to override
  header-based dispatch.

- index:

  Integer volume index. Used as the single-volume selector for
  [`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
  (when dispatching to it) and, when `indices` is `NULL`, forwarded to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
  as `indices` so you can pull out a subset of time points while still
  receiving a `NeuroVec`.

- indices:

  Optional integer vector of sub-volume indices forwarded to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md).
  Takes precedence over `index`.

- mask:

  Optional spatial mask passed through to the vector or hyper-vector
  readers. Ignored by `read_vol`.

- mode:

  IO mode forwarded to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md);
  see that function for details. Ignored for 3D and 5D dispatch.

## Value

The return type depends on dispatch:

- **3D dispatch** (single effectively-3D file, or `type = "vol"`): a
  [`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md).

- **4D dispatch** (single 4D file, or `type = "vec"` with one file): a
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  (concrete subclass depends on `mode`).

- **Multi-file dispatch** (`length(file_name) > 1`, or `type = "vec"`
  with multiple files): a
  [`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md),
  which extends `NeuroVec` but stores each file as a distinct segment
  rather than concatenating into a single contiguous 4D array. Mixed
  3D/4D inputs are allowed; each 3D file contributes one time point to
  the sequence.

- **5D dispatch** (single 5D file, or `type = "hyper"`): a
  [`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md).

## Details

Auto-dispatch (`type = "auto"`) uses the following rules:

- `length(file_name) > 1`: always routed to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md),
  so the result is a
  [`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md)
  regardless of whether the individual files are 3D, 4D, or a mix. See
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
  for the return-type details.

- Single file with a 5th dimension `> 1`: routed to
  [`read_hyper_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_hyper_vec.md),
  returning a
  [`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md).

- Single file with a 4th dimension `> 1`: routed to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md),
  returning a
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md).
  If `index` is supplied (and `indices` is not), it is forwarded as
  `indices` so you can pull out a subset while still getting a
  `NeuroVec` back.

- Single effectively-3D file (either truly 3D or 4D with `dim[4] == 1`):
  routed to
  [`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md),
  returning a
  [`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md).

Explicit `type` values bypass header inspection:

- `type = "vol"`: requires a single file; always returns a
  [`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md).
  The `index` argument picks the sub-volume when the file is 4D.

- `type = "vec"`: forwards to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md).
  Multi-file input yields a
  [`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md).
  A single 3D file is promoted to a `NeuroVec` with `dim[4] == 1`.

- `type = "hyper"`: requires a single file; returns a
  [`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md).

## See also

[`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md),
[`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md),
[`read_hyper_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_hyper_vec.md)

## Examples

``` r
# 3D file -> DenseNeuroVol
vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))

# 4D file -> NeuroVec
vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Multiple files (any mix of 3D and 4D) -> NeuroVecSeq
fn <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
seq_vec <- read_image(c(fn, fn))
class(seq_vec)  # "NeuroVecSeq"
#> [1] "NeuroVecSeq"
#> attr(,"package")
#> [1] "neuroim2"
```
