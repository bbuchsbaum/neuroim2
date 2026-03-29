# Leave-one-run-out smoothing helper

Leave-one-run-out smoothing helper

## Usage

``` r
cgb_smooth_loro(runs, graphs, passes = 1L, lambda = 1)
```

## Arguments

- runs:

  List of
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  objects (one per run).

- graphs:

  List of graphs returned by `cgb_make_graph(..., leave_one_out=TRUE)`.

- passes, lambda:

  See
  [`cgb_smooth`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_smooth.md).

## Value

A list of smoothed `NeuroVec` objects, one per run.
