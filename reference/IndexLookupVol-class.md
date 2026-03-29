# IndexLookupVol Class

A three-dimensional brain image class that serves as a map between 1D
grid indices and a table of values. This class is primarily used in
conjunction with the
[`SparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
class to efficiently represent and access sparse neuroimaging data.

The `IndexLookupVol` class provides efficient indexing and coordinate
lookup functionality for 3D neuroimaging data. It maintains a mapping
between linear indices and 3D coordinates, optimizing memory usage and
access speed for sparse volumes.

Creates an `IndexLookupVol` object, which provides efficient
bidirectional mapping between linear indices and 3D coordinates in a
neuroimaging volume. This is particularly useful for working with masked
or sparse brain volumes.

## Usage

``` r
IndexLookupVol(space, indices)
```

## Arguments

- space:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the 3D space dimensions, spacing, and orientation.

- indices:

  An integer vector containing the linear indices of the voxels to
  include in the lookup volume. These should be 1-based indices within
  the range of the space.

## Value

An object of class `IndexLookupVol` containing:

- A mapping between linear indices and sparse positions

- The original space information

- The subset of included voxel indices

## Details

The IndexLookupVol class extends
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
and provides a mechanism for efficient lookup and mapping of sparse 3D
neuroimaging data. It stores only the indices of non-zero voxels and
their corresponding mappings, allowing for memory-efficient
representation of large, sparse brain images.

Create an IndexLookupVol Object

## Slots

- `space`:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object representing the 3D space of the brain image.

- `indices`:

  An integer vector containing the 1D indices of the non-zero voxels in
  the grid.

- `map`:

  An integer vector containing the mapping between the 1D indices and
  the table of values.

## Methods

This class inherits methods from
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md).
Additional methods specific to index lookup and mapping operations may
be available.

## Implementation Details

The class uses an integer mapping array for O(1) lookups between linear
indices and their corresponding positions in the sparse representation.

## See also

[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for the primary class that utilizes IndexLookupVol.
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
for the base volumetric image class.

`IndexLookupVol` for creating instances of this class

[`coords`](https://bbuchsbaum.github.io/neuroim2/reference/coords.md)
for coordinate lookup,
[`lookup`](https://bbuchsbaum.github.io/neuroim2/reference/lookup.md)
for index mapping,
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for space representation

## Examples

``` r
# Create a NeuroSpace object
space <- NeuroSpace(dim = c(2L, 2L, 2L), origin = c(0, 0, 0), spacing = c(1, 1, 1))

# Create a 3D mask
mask <- array(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE), dim = c(2, 2, 2))

# Create indices and map for the IndexLookupVol
indices <- which(mask)
map <- seq_along(indices)

# Create an IndexLookupVol object
ilv <- IndexLookupVol(space = space, indices = as.integer(indices))

# Access the indices
print(ilv@indices)
#> [1] 1 3 6 8

# Access the map
print(ilv@map)
#> [1] 1 0 2 0 0 3 0 4


# Create a 64x64x64 space
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))

# Create a lookup volume with random indices
indices <- sample(1:262144, 10000)  # Select 10000 random voxels
ilv <- IndexLookupVol(space, indices)

# Look up coordinates for specific indices
coords <- coords(ilv, indices[1:10])

```
