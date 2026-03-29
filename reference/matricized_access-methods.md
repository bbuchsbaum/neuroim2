# Extract values from a 4D tensor using a matrix of time-space indices.

This function efficiently extracts values from a 4D tensor (typically
neuroimaging data) using a matrix of indices where each row contains a
time index in column 1 and a spatial index in column 2. The spatial
index refers to the position in the flattened spatial dimensions
(x,y,z). This is primarily used internally by the
[`series()`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md)
method to efficiently access time series data for specific voxels.

## Usage

``` r
matricized_access(x, i, ...)

# S4 method for class 'SparseNeuroVec,matrix'
matricized_access(x, i)

# S4 method for class 'SparseNeuroVec,integer'
matricized_access(x, i)

# S4 method for class 'SparseNeuroVec,numeric'
matricized_access(x, i)

# S4 method for class 'BigNeuroVec,matrix'
matricized_access(x, i)

# S4 method for class 'BigNeuroVec,integer'
matricized_access(x, i)

# S4 method for class 'BigNeuroVec,numeric'
matricized_access(x, i)
```

## Arguments

- x:

  a data source, typically a `SparseNeuroVec` object containing 4D
  neuroimaging data

- i:

  Either:

  - A matrix with 2 columns: \[time_index, space_index\] specifying
    which values to extract

  - A numeric vector of spatial indices to extract all timepoints for
    those locations

- ...:

  additional arguments to be passed to methods.

## Value

When `i` is a matrix, returns a numeric vector of values at the
specified time-space coordinates. When `i` is a vector, returns a matrix
where each column contains the full time series for each spatial index.

## Examples

``` r
# Create a sparse 4D neuroimaging vector
bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)

# Extract specific timepoint-voxel pairs
# Get value at timepoint 1, voxel 1 and timepoint 2, voxel 2
idx_mat <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
vals <- matricized_access(svec, idx_mat)

# Get full time series for voxels 1 and 2
ts_mat <- matricized_access(svec, c(1,2))
# Each column in ts_mat contains the full time series for that voxel
```
