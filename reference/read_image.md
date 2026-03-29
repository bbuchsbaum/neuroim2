# read_image

Convenience wrapper that inspects the file metadata and dispatches to
[`read_vol`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
for 3D data,
[`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
for 4D data, or
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

  Character vector of file paths.

- type:

  One of `"auto"`, `"vol"`, `"vec"`, or `"hyper"` to override dispatch.

- index:

  Volume index to use when returning a
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  or when you want to load a subset of volumes while still returning a
  `NeuroVec`.

- indices:

  Optional vector of indices passed through to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md).

- mask:

  Optional spatial mask passed through to vector/hyper-vector readers.

- mode:

  IO mode forwarded to
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md).

## Value

A
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
when the input is effectively 3D (or when `type = "vol"`), a
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)/[`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md)
for 4D input, or a
[`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md)
for 5D input.

## Examples

``` r
vol <- read_image(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
vec <- read_image(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
```
