# SparseNeuroVecSource Class

A class used to produce a
[`SparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
instance. It encapsulates the necessary information to create a sparse
representation of a 4D neuroimaging dataset.

## Details

SparseNeuroVecSource acts as a factory for SparseNeuroVec objects. It
holds the spatial mask that determines which voxels will be included in
the sparse representation. This class is typically used in data loading
or preprocessing pipelines where the sparse structure of the data is
known or determined before the full dataset is loaded.

## Slots

- `mask`:

  An object of class
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  representing the subset of voxels that will be stored in memory. This
  mask defines the sparse structure of the resulting SparseNeuroVec.

## Inheritance

`SparseNeuroVecSource` inherits from:

- [`NeuroVecSource`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSource-class.md):
  Base class for NeuroVec source objects

## See also

[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for the resulting sparse 4D neuroimaging data class.
[`LogicalNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
for the mask representation.

## Examples

``` r
# Create a simple mask
mask_data <- array(runif(64*64*32) > 0.7, dim = c(64, 64, 32))
mask <- LogicalNeuroVol(mask_data, space = NeuroSpace(dim = c(64, 64, 32)))

# Create a SparseNeuroVecSource
sparse_source <- new("SparseNeuroVecSource", mask = mask)

```
