# Load image data from a NeuroVecSource object

This function loads the image data from a NeuroVecSource object,
handling various dimensionalities and applying any necessary
transformations.

## Usage

``` r
# S4 method for class 'MappedNeuroVecSource'
load_data(x)

# S4 method for class 'NeuroVecSource'
load_data(x)

# S4 method for class 'NeuroVolSource'
load_data(x)

# S4 method for class 'SparseNeuroVecSource'
load_data(x)
```

## Arguments

- x:

  The NeuroVecSource object containing the image metadata and file
  information.

## Value

a DenseNeuroVec object

## Details

This method performs the following steps: 1. Validates the
dimensionality of the metadata. 2. Reads the image data using RNifti. 3.
Handles 5D arrays by dropping the 4th dimension if it has length 1. 4.
Applies slope scaling if present in the metadata. 5. Constructs a
NeuroSpace object with appropriate dimensions and spatial information.
6. Creates and returns a DenseNeuroVec object, handling both 3D and 4D
input arrays.

## Note

This method currently only supports NIfTI file format through RNifti.

## See also

[`NeuroVecSource`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSource.md),
[`DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md),
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace.md)
