# Extract Mask from Neuroimaging Object

Generic function to extract or generate a mask from neuroimaging
objects. For sparse objects with a `@mask` slot, returns the stored
mask. For dense objects, returns a filled mask (all TRUE values)
indicating all voxels contain valid data.

## Usage

``` r
mask(x)

# S4 method for class 'ClusteredNeuroVol'
mask(x)

# S4 method for class 'FileBackedNeuroVec'
mask(x)

# S4 method for class 'MappedNeuroVec'
mask(x)

# S4 method for class 'NeuroHyperVec'
mask(x)

# S4 method for class 'NeuroSlice'
mask(x)

# S4 method for class 'DenseNeuroVec'
mask(x)

# S4 method for class 'DenseNeuroVol'
mask(x)

# S4 method for class 'LogicalNeuroVol'
mask(x)

# S4 method for class 'AbstractSparseNeuroVec'
mask(x)

# S4 method for class 'SparseNeuroVecSource'
mask(x)
```

## Arguments

- x:

  A neuroimaging object (NeuroVol, NeuroVec, or derived classes)

## Value

A
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
object representing the mask

## Details

The behavior depends on the class of the input object:

- For sparse objects (SparseNeuroVec, ClusteredNeuroVol, etc.): Returns
  the stored `@mask` slot

- For dense objects (DenseNeuroVol, DenseNeuroVec, etc.): Returns a
  LogicalNeuroVol with all TRUE values

- For ROI objects: Not implemented (use
  [`coords()`](https://bbuchsbaum.github.io/neuroim2/reference/coords.md)
  instead)

## Examples

``` r
# Create a dense volume
vol <- NeuroVol(array(rnorm(64^3), c(64,64,64)), NeuroSpace(c(64,64,64)))
m <- mask(vol)  # Returns all TRUE mask

# Create a sparse vector with explicit mask
mask_array <- array(runif(64^3) > 0.5, c(64,64,64))
mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(c(64,64,64)))
# Data must be a matrix (time x masked voxels)
sparse_data <- matrix(rnorm(sum(mask_array) * 10), nrow = 10, ncol = sum(mask_array))
svec <- SparseNeuroVec(sparse_data, NeuroSpace(c(64,64,64,10)), mask_vol)
m2 <- mask(svec)  # Returns the stored mask
```
