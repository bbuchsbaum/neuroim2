# Reorient Image to Canonical (RAS+) Orientation

Reorients a neuroimaging volume or vector to the canonical RAS+
(Right-Anterior-Superior) orientation by permuting and flipping axes.
This is equivalent to nibabel's `as_closest_canonical()`.

## Usage

``` r
as_canonical(x, target = c("R", "A", "S"))
```

## Arguments

- x:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  or
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  object.

- target:

  Character vector of length 3 giving the desired orientation. Default
  is `c("R", "A", "S")` (RAS+).

## Value

A reoriented object of the same class as `x`.

## Details

The function works by computing the orientation transform from the
current axis codes to the target codes, then applying the necessary axis
permutations and flips. The affine matrix is updated to reflect the new
orientation while preserving world-coordinate mapping.

## See also

[`axcodes`](https://bbuchsbaum.github.io/neuroim2/reference/axcodes-methods.md),
[`reorient`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)

## Examples

``` r
sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
vol <- DenseNeuroVol(array(rnorm(1000), c(10,10,10)), sp)
ras_vol <- as_canonical(vol)
axcodes(ras_vol)  # "R" "A" "S"
#> [1] "R" "A" "S"
```
