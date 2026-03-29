# MappedNeuroVecSource Class

A class used to produce a
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
instance. It encapsulates the necessary information to create a
memory-mapped representation of a 4D neuroimaging dataset.

Creates a `MappedNeuroVecSource` object that manages the memory mapping
between a neuroimaging file and memory space. This is typically used
internally by
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
but can be created directly for custom access patterns.

## Usage

``` r
MappedNeuroVecSource(file_name)
```

## Arguments

- file_name:

  Character string specifying the path to the neuroimaging file.
  Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).

## Value

A new `MappedNeuroVecSource` object containing:

- Meta information about the dataset

- File format details

- Dimensional information

## Details

MappedNeuroVecSource acts as a factory for MappedNeuroVec objects. While
it doesn't have any additional slots beyond its parent class, it
specifies the intent to create a memory-mapped representation of the
neuroimaging data. This class is typically used in data loading
pipelines where large datasets need to be accessed efficiently without
loading the entire dataset into memory.

Create a Memory-Mapped Source for Neuroimaging Data

The function performs several important checks:

- Validates file existence and permissions

- Reads and validates header information

- Ensures proper dimensionality (\>= 3D)

- Verifies file format compatibility

## Inheritance

`MappedNeuroVecSource` inherits from:

- [`NeuroVecSource`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSource-class.md):
  Base class for NeuroVec source objects

## See also

[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
for the main user interface,
[`read_header`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md)
for header reading details

## Examples

``` r
# Create a MappedNeuroVecSource
mapped_source <- new("MappedNeuroVecSource")

# \donttest{
# Create source from NIFTI file
source <- MappedNeuroVecSource(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Check dimensions
dim(source@meta_info)
#> [1] 64 64 25  4

# View header information
str(source@meta_info)
#> Formal class 'NIFTIMetaInfo' [package "neuroim2"] with 17 slots
#>   ..@ nifti_header     : list()
#>   ..@ header_file      : chr "/home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii"
#>   ..@ data_file        : chr "/home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii"
#>   ..@ descriptor       :Formal class 'NIFTIFormat' [package "neuroim2"] with 5 slots
#>   .. .. ..@ file_format     : chr "NIFTI"
#>   .. .. ..@ header_encoding : chr "raw"
#>   .. .. ..@ header_extension: chr "nii"
#>   .. .. ..@ data_encoding   : chr "raw"
#>   .. .. ..@ data_extension  : chr "nii"
#>   ..@ endian           : chr "little"
#>   ..@ data_offset      : num 352
#>   ..@ bytes_per_element: int 4
#>   ..@ intercept        : num 0
#>   ..@ slope            : num 1
#>   ..@ header           :List of 42
#>   .. ..$ file_type     : chr "NIfTI"
#>   .. ..$ encoding      : chr "binary"
#>   .. ..$ version       : chr "1"
#>   .. ..$ file_name     : chr "/home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii"
#>   .. ..$ endian        : chr "little"
#>   .. ..$ diminfo       : int 0
#>   .. ..$ dimensions    : num [1:8] 4 64 64 25 4 1 1 1
#>   .. ..$ num_dimensions: num 4
#>   .. ..$ intent1       : num 0
#>   .. ..$ intent2       : num 0
#>   .. ..$ intent3       : num 0
#>   .. ..$ intent_code   : int 0
#>   .. ..$ datatype      : int 16
#>   .. ..$ data_storage  : chr "FLOAT"
#>   .. ..$ bitpix        : int 32
#>   .. ..$ slice_start   : int 0
#>   .. ..$ pixdim        : num [1:8] -1 3.5 3.5 3.7 0 ...
#>   .. ..$ qfac          : num -1
#>   .. ..$ vox_offset    : num 352
#>   .. ..$ scl_slope     : num 1
#>   .. ..$ scl_intercept : num 0
#>   .. ..$ slice_end     : int 0
#>   .. ..$ slice_code    : int 0
#>   .. ..$ xyzt_units    : int 2
#>   .. ..$ cal_max       : num 0
#>   .. ..$ cal_min       : num 0
#>   .. ..$ slice_duration: num 0
#>   .. ..$ toffset       : num 0
#>   .. ..$ glmax         : int 0
#>   .. ..$ glmin         : int 0
#>   .. ..$ description   : int [1:80] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ auxfile       : int [1:24] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ qform_code    : int 1
#>   .. ..$ sform_code    : int 1
#>   .. ..$ quaternion    : num [1:3] 0 1 0
#>   .. ..$ qoffset       : num [1:3] 112 -108 -46.2
#>   .. ..$ qform         : num [1:4, 1:4] -3.5 0 0 0 0 3.5 0 0 0 0 ...
#>   .. ..$ sform         : num [1:4, 1:4] -3.5 0 0 0 0 3.5 0 0 0 0 ...
#>   .. ..$ intent_name   : chr [1:16] "" "" "" "" ...
#>   .. ..$ magic         : chr "n+1"
#>   .. ..$ onefile       : logi TRUE
#>   .. ..$ extensions    :Formal class 'NiftiExtensionList' [package "neuroim2"] with 1 slot
#>   .. .. .. ..@ .Data: list()
#>   ..@ data_type        : chr "FLOAT"
#>   ..@ dims             : num [1:4] 64 64 25 4
#>   ..@ spatial_axes     :Formal class 'AxisSet3D' [package "neuroim2"] with 4 slots
#>   .. .. ..@ k   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. ..@ axis     : chr "Inferior-to-Superior"
#>   .. .. .. .. ..@ direction: num [1:3] 0 0 1
#>   .. .. ..@ j   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. ..@ axis     : chr "Posterior-to-Anterior"
#>   .. .. .. .. ..@ direction: num [1:3] 0 1 0
#>   .. .. ..@ i   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. ..@ axis     : chr "Right-to-Left"
#>   .. .. .. .. ..@ direction: num [1:3] -1 0 0
#>   .. .. ..@ ndim: int 3
#>   ..@ additional_axes  :Formal class 'AxisSet' [package "neuroim2"] with 1 slot
#>   .. .. ..@ ndim: int 0
#>   ..@ spacing          : num [1:3] 3.5 3.5 3.7
#>   ..@ origin           : num [1:3] 112 -108 -46.2
#>   ..@ label            : chr "global_mask_v4"
# }
```
