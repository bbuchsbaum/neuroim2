# Convert to dense representation

Convert to dense representation

## Usage

``` r
as.dense(x)
```

## Arguments

- x:

  the object to densify

## Value

A dense representation of the input object.

## Examples

``` r
# Create a sparse representation
space <- NeuroSpace(c(10,10,10,4), c(1,1,1))
mask <- array(runif(10*10*10) > 0.8, c(10,10,10))  # ~20% of voxels active
data <- matrix(rnorm(sum(mask) * 4), 4, sum(mask))  # Random data for active voxels
sparse_vec <- SparseNeuroVec(data, space, mask)

# Convert to dense representation
dense_vec <- as.dense(sparse_vec)
# The dense representation has the same dimensions but stores all voxels
identical(dim(sparse_vec), dim(dense_vec))
#> [1] TRUE
```
