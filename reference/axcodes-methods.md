# Get Orientation Axis Codes

Returns a character vector of anatomical axis direction labels for a
neuroimaging object or space. For example, `c("R", "A", "S")` for a
standard RAS-oriented image.

## Usage

``` r
axcodes(x)

# S4 method for class 'NeuroSpace'
axcodes(x)

# S4 method for class 'NeuroObj'
axcodes(x)

# S4 method for class 'matrix'
axcodes(x)
```

## Arguments

- x:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md),
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md),
  or a 4x4 affine matrix.

## Value

A character vector of length 3 with axis direction codes.

## Examples

``` r
sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
axcodes(sp)
#> [1] "R" "A" "S"

vol <- DenseNeuroVol(array(0, c(10,10,10)), sp)
axcodes(vol)
#> [1] "R" "A" "S"
```
