# read_vol_list

This function loads a list of image volumes and returns a NeuroVec
object.

## Usage

``` r
read_vol_list(file_names, mask = NULL)
```

## Arguments

- file_names:

  A list of file names to load.

- mask:

  An optional mask defining the subset of voxels to load.

## Value

An instance of the
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
class.
