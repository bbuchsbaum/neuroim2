# Apply a precomputed CGB graph to volumetric data

Apply a precomputed CGB graph to volumetric data

## Usage

``` r
cgb_smooth(x, graph, passes = 1L, lambda = 1)
```

## Arguments

- x:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  (4D) or
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  (3D).

- graph:

  Graph list returned by
  [`cgb_make_graph`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_make_graph.md).

- passes:

  Number of smoothing passes (\>=1). Each pass multiplies by `W`; if
  `lambda < 1` a simple diffusion blend `(1 - lambda)I + lambda W` is
  applied per pass.

- lambda:

  Blend factor in \\\[0,1\]\\ controlling diffusion strength.

## Value

Smoothed object of the same class as `x`.
