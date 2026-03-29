# Convert to a LogicalNeuroVol

Convert to a LogicalNeuroVol

## Usage

``` r
as.mask(x, indices)
```

## Arguments

- x:

  the object to binarize

- indices:

  the indices to set to TRUE

## Value

A `LogicalNeuroVol` object with `TRUE` values at the specified
`indices`.

## Examples

``` r
# Create a simple 3D volume with random values
space <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
vol <- NeuroVol(array(runif(1000), c(10,10,10)), space)

# Create a mask by thresholding (values > 0.5 become TRUE)
mask1 <- as.mask(vol > 0.5)

# Create a mask by specifying indices
indices <- which(vol > 0.8)  # get indices of high values
mask2 <- as.mask(vol, indices)

# Both masks are LogicalNeuroVol objects
identical(class(mask1), class(mask2))
#> [1] TRUE
```
