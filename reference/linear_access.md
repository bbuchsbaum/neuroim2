# Extract values from an array-like object using linear indexing.

This function extracts the values of the elements in an array-like
object using linear indexing. Linear indexing is a way of indexing an
array by a single index that is computed from multiple indices using a
formula.

## Usage

``` r
linear_access(x, i, ...)
```

## Arguments

- x:

  a data source.

- i:

  a vector of indices.

- ...:

  additional arguments to be passed to methods.

## Value

A `vector` containing the values at the specified linear indices of `x`.

## Examples

``` r
# Create a sparse neuroimaging vector
bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)

# Extract values using linear indices
# Get values from first timepoint at voxels 1,2,3
indices <- c(1,2,3)
vals <- linear_access(svec, indices)

# Get values from multiple timepoints and voxels
# First voxel at timepoint 1, second voxel at timepoint 2
indices <- c(1, 1000 + 2) # 1000 = prod(10,10,10)
vals <- linear_access(svec, indices)
```
