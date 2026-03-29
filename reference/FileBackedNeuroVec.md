# Create a File-Backed Neuroimaging Vector

Constructs a
[`FileBackedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/FileBackedNeuroVec-class.md)
instance, which represents a file-backed neuroimaging vector object.
This constructor provides memory-efficient access to large neuroimaging
datasets by keeping the data on disk until needed.

## Usage

``` r
FileBackedNeuroVec(file_name, label = basename(file_name))
```

## Arguments

- file_name:

  A character string specifying the path to the neuroimaging file.
  Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).

- label:

  Optional character string providing a label for the vector

## Value

A new instance of class
[`FileBackedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/FileBackedNeuroVec-class.md).

## Details

Create a FileBackedNeuroVec Object

The function performs the following operations:

- Reads the header information from the specified file

- Validates the dimensionality (must be 4D data)

- Creates a
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object with appropriate metadata

- Initializes the file-backed vector with minimal memory footprint

## See also

[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for spatial metadata management,
[`read_header`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md)
for header information extraction,
[`sub_vector`](https://bbuchsbaum.github.io/neuroim2/reference/sub_vector-methods.md)
for data access methods

## Examples

``` r
# Create a file-backed vector from a NIFTI file
fbvec <- FileBackedNeuroVec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Access specific volumes without loading entire dataset
first_vol <- sub_vector(fbvec, 1)

```
