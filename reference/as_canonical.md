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

Reorienting to a target that differs only in which anatomical direction
each axis runs is a relabelling of the grid, not a resampling: the
voxels are the same voxels, permuted and flipped. This function does
exactly that, so it is exact, works on images of any size, and costs a
single [`aperm()`](https://rdrr.io/r/base/aperm.html).

It used to build the target space and hand it to
[`resample`](https://bbuchsbaum.github.io/neuroim2/reference/resample-methods.md),
which had three consequences: values were interpolated where they should
merely have moved, the registration backend refused images smaller than
four voxels on any axis, and – because the target space kept the
source's dimensions – an axis-permuting reorientation wrote into a grid
that did not contain the data and silently lost about a quarter of it.

## See also

[`axcodes`](https://bbuchsbaum.github.io/neuroim2/reference/axcodes-methods.md),
[`reorient`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md),
[`apply_orientation`](https://bbuchsbaum.github.io/neuroim2/reference/orientation_utils.md)

## Examples

``` r
sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
vol <- DenseNeuroVol(array(rnorm(1000), c(10,10,10)), sp)
ras_vol <- as_canonical(vol)
axcodes(ras_vol)  # "R" "A" "S"
#> [1] "R" "A" "S"
```
