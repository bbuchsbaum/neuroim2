# BigNeuroVec Class

A class representing a sparse four-dimensional brain image backed by a
disk-based big matrix. BigNeuroVec objects are designed for efficient
handling of large-scale brain imaging data that exceeds available
memory.

## Details

BigNeuroVec leverages file-backed storage to manage large 4D
neuroimaging datasets that would typically exceed available RAM. It
combines the sparse representation framework of
[`AbstractSparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md)
with the disk-based storage capabilities of `FBM`, allowing for
out-of-core computations on massive datasets.

## Slots

- `data`:

  An instance of class `FBM` from the `bigstatsr` package, containing
  time-series data. The FBM (File-Backed Big Matrix) is a matrix-like
  structure stored on disk, enabling efficient handling of large-scale
  data.

## Inheritance

`BigNeuroVec` inherits from:

- [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md):
  Base class for 4D brain images

- [`AbstractSparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md):
  Provides sparse representation framework

- `ArrayLike4D`: Interface for 4D array-like operations

## See also

[`AbstractSparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md)
for the parent sparse representation class.
[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the base 4D brain image class.
[`FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
for details on File-Backed Big Matrix objects.
