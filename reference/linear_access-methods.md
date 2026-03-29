# Linear Access Method for FileBackedNeuroVec

Internal method providing linear access to memory-mapped data.

Provides linear access to the data across all vectors in the sequence.

## Usage

``` r
# S4 method for class 'DenseNeuroVol,numeric'
linear_access(x, i)

# S4 method for class 'DenseNeuroVec,numeric'
linear_access(x, i)

# S4 method for class 'DenseNeuroVol,integer'
linear_access(x, i)

# S4 method for class 'DenseNeuroVec,integer'
linear_access(x, i)

# S4 method for class 'FileBackedNeuroVec,numeric'
linear_access(x, i)

# S4 method for class 'MappedNeuroVec,numeric'
linear_access(x, i)

# S4 method for class 'NeuroHyperVec,ANY'
linear_access(x, i, ...)

# S4 method for class 'NeuroVecSeq,numeric'
linear_access(x, i)

# S4 method for class 'SparseNeuroVol,numeric'
linear_access(x, i)

# S4 method for class 'AbstractSparseNeuroVec,numeric'
linear_access(x, i)
```

## Arguments

- x:

  A NeuroVecSeq object

- i:

  Numeric vector of indices for linear access

- ...:

  Additional arguments (not used)

## Value

Numeric vector of accessed values

## Examples

``` r
# \donttest{
# Create a small NeuroVec and save it
nvec <- NeuroVec(matrix(1:32, 8, 4), NeuroSpace(c(2,2,2,4)))
tmp <- tempfile(fileext = ".nii")
write_vec(nvec, tmp)

# Load as FileBackedNeuroVec and access values
fbvec <- FileBackedNeuroVec(tmp)
values <- linear_access(fbvec, 1:10)

# Clean up
unlink(tmp)
# }
```
