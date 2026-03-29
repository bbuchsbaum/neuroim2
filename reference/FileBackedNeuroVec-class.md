# FileBackedNeuroVec Class

A class representing a four-dimensional brain image that uses on-demand
loading through memory-mapped file access. This approach enables
efficient handling of large-scale brain imaging data by loading only the
required portions of the data into memory when needed.

The `FileBackedNeuroVec` class represents a memory-efficient vector of
neuroimaging data that is stored on disk rather than in memory. This is
particularly useful for large datasets where memory constraints are a
concern.

## Details

FileBackedNeuroVec objects provide a memory-efficient solution for
working with large 4D neuroimaging datasets. By utilizing memory-mapped
file access, this class allows users to work with datasets that exceed
available RAM, only loading the necessary data segments into memory as
they are accessed.

## Slots

- `meta`:

  An instance of class
  [`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)
  containing file metadata such as file path, format, and other
  associated information.

## Inheritance

`FileBackedNeuroVec` inherits from:

- [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md):
  Base class for 4D brain images

- `ArrayLike4D`: Interface for 4D array-like operations

## Memory Management

Data is read from disk on-demand, reducing memory usage compared to
in-memory storage. The trade-off is slightly slower access times due to
disk I/O operations.

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the base 4D brain image class.
[`FileMetaInfo-class`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)
for details on file metadata representation.

[`FileBackedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/FileBackedNeuroVec.md)
for creating instances of this class

## Examples

``` r
# Load example 4D image file included with package
file_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
fbvec <- FileBackedNeuroVec(file_path)

# Get dimensions of the image
dim(fbvec)
#> [1] 64 64 25  4

# Extract first volume
vol1 <- sub_vector(fbvec, 1)

# Extract multiple volumes
vols <- sub_vector(fbvec, 1:2)
```
