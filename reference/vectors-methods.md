# Extract an ordered list of 1D vectors.

This function extracts an ordered list of 1D vectors from an object that
supplies vector data. The `subset` argument specifies the subset of
vectors to extract, and can be a vector of indices or a logical vector.
The return value is a list containing the extracted vectors in the same
order as the specified indices.

## Usage

``` r
vectors(x, subset, ...)

# S4 method for class 'NeuroVec,missing'
vectors(x)

# S4 method for class 'DenseNeuroVec,missing'
vectors(x)

# S4 method for class 'NeuroVec,numeric'
vectors(x, subset)

# S4 method for class 'NeuroVec,logical'
vectors(x, subset)

# S4 method for class 'NeuroVecSeq,missing'
vectors(x)

# S4 method for class 'NeuroVecSeq,numeric'
vectors(x, subset)

# S4 method for class 'NeuroVecSeq,logical'
vectors(x, subset)

# S4 method for class 'ROIVec,missing'
vectors(x)

# S4 method for class 'matrix,missing'
vectors(x)

# S4 method for class 'ROIVec,integer'
vectors(x, subset)

# S4 method for class 'matrix,integer'
vectors(x, subset)

# S4 method for class 'matrix,numeric'
vectors(x, subset)

# S4 method for class 'ROIVec,numeric'
vectors(x, subset)

# S4 method for class 'ROIVec,logical'
vectors(x, subset)

# S4 method for class 'SparseNeuroVec,missing'
vectors(x, nonzero = FALSE)
```

## Arguments

- x:

  the object that supplies the vector data.

- subset:

  the subset of vectors to extract.

- ...:

  additional arguments to be passed to methods.

- nonzero:

  only include nonzero vectors in output list

## Value

A `list` containing the extracted vectors from `x` in the same order as
`subset`.

A deflist object where each element is a function that returns the time
series for a voxel. The length of the deflist equals the total number of
voxels.

## Examples

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vec <- read_vec(file_name)
v <- vectors(vec)
mean(v[[1]])
#> [1] 0
```
