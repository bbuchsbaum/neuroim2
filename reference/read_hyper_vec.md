# Read a 5D image as a NeuroHyperVec

Read a 5D image as a NeuroHyperVec

## Usage

``` r
read_hyper_vec(file_name, mask = NULL)
```

## Arguments

- file_name:

  Path to a single NIfTI file.

- mask:

  Optional spatial mask (logical array/vector, `NeuroVol`, or
  `LogicalNeuroVol`). When provided, only masked voxels are stored.

## Value

A
[`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md).
