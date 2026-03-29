# AbstractSparseNeuroVec Class

An abstract base class for sparse four-dimensional brain image
representations. This class provides the foundation for efficient
storage and manipulation of large, sparse neuroimaging data.

## Details

The AbstractSparseNeuroVec class serves as a template for implementing
various sparse representations of 4D brain images. It combines the
spatial properties of
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
with the efficiency of sparse data structures.

## Slots

- `mask`:

  An object of class
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  defining the sparse domain of the brain image. This mask indicates
  which voxels contain non-zero data.

- `map`:

  An object of class
  [`IndexLookupVol`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
  used to map between spatial coordinates and index/row coordinates in
  the sparse representation.

## Subclasses

Concrete implementations of this abstract class should provide specific
data storage mechanisms and methods for efficient access and
manipulation of sparse 4D brain image data.

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the parent class.
[`LogicalNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
for the mask representation.
[`IndexLookupVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
for the spatial-to-index mapping.
