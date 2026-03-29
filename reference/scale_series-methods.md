# Generic functions to scale (center and/or normalize by standard deviation) each series of a 4D image That is, if the 4th dimension is 'time' each series is a 1D time series.

Generic functions to scale (center and/or normalize by standard
deviation) each series of a 4D image That is, if the 4th dimension is
'time' each series is a 1D time series.

## Usage

``` r
scale_series(x, center, scale)

# S4 method for class 'NeuroVec,logical,missing'
scale_series(x, center, scale)

# S4 method for class 'DenseNeuroVec,logical,logical'
scale_series(x, center, scale)

# S4 method for class 'SparseNeuroVec,logical,logical'
scale_series(x, center, scale)

# S4 method for class 'NeuroVec,logical,logical'
scale_series(x, center, scale)

# S4 method for class 'NeuroVec,missing,logical'
scale_series(x, center, scale)

# S4 method for class 'NeuroVec,missing,missing'
scale_series(x, center, scale)
```

## Arguments

- x:

  a four dimensional image

- center:

  a `logical` value indicating whether series should be centered

- scale:

  a `logical` value indicating whether series should be divided by
  standard deviation

## Value

An object of the same class as `x`, with each time series centered
and/or scaled.

## Examples

``` r
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
res <- scale_series(bvec, TRUE, TRUE)
```
