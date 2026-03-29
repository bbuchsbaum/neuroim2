# Create Neuroimaging Metadata Object

Creates a MetaInfo object containing essential metadata for neuroimaging
data, including dimensions, spacing, orientation, and data type
information.

## Usage

``` r
MetaInfo(
  Dim,
  spacing,
  origin = rep(0, length(spacing)),
  data_type = "FLOAT",
  label = "",
  spatial_axes = OrientationList3D$AXIAL_LPI,
  additional_axes = NullAxis
)
```

## Arguments

- Dim:

  Integer vector. Image dimensions (e.g., c(64, 64, 32) for 3D).

- spacing:

  Numeric vector. Voxel dimensions in mm.

- origin:

  Numeric vector. Coordinate origin. Default is zero vector.

- data_type:

  Character. Data type (e.g., "FLOAT", "SHORT"). Default is "FLOAT".

- label:

  Character. Image label(s). Default is "".

- spatial_axes:

  Object. Spatial orientation. Default is OrientationList3D\$AXIAL_LPI.

- additional_axes:

  Object. Non-spatial axes. Default is NullAxis.

## Value

A MetaInfo object

## Details

Create MetaInfo Object

The MetaInfo object is fundamental for:

- Spatial interpretation of image data

- Data type handling and conversion

- Memory allocation and mapping

- File I/O operations

Input validation ensures:

- Dimensions are positive integers

- Spacing values are positive

- Origin coordinates are finite

- Data type is supported

## See also

[`NIFTIMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md),
[`AFNIMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)

## Examples

``` r
# Create metadata for 3D structural MRI
meta <- MetaInfo(
  Dim = c(256, 256, 180),
  spacing = c(1, 1, 1),
  data_type = "FLOAT",
  label = "T1w"
)

# Get image dimensions
dim(meta)
#> NULL

# Get transformation matrix
trans(meta)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    0    1
```
