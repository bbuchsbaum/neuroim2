# MetaInfo Class

This class encapsulates meta information for neuroimaging data types,
including spatial and temporal characteristics, data type, and labeling.

## Details

The MetaInfo class provides a structured way to store and access
essential metadata for neuroimaging data. This includes information
about the data type, spatial and temporal dimensions, voxel spacing, and
coordinate system origin.

## Slots

- `data_type`:

  A `character` string specifying the data type code (e.g., "FLOAT",
  "INT").

- `dims`:

  A `numeric` vector representing image dimensions.

- `spatial_axes`:

  An
  [`AxisSet3D`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet3D-class.md)
  object representing image axes for spatial dimensions (x, y, z).

- `additional_axes`:

  An
  [`AxisSet`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)
  object representing axes for dimensions beyond spatial (e.g., time,
  color band, direction).

- `spacing`:

  A `numeric` vector representing voxel dimensions in real-world units.

- `origin`:

  A `numeric` vector representing the coordinate origin.

- `label`:

  A `character` vector containing name(s) of images or data series.

## See also

[`FileMetaInfo-class`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md),
[`AxisSet3D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet3D-class.md),
[`AxisSet-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)

## Examples

``` r
# Create a MetaInfo object
meta_info <- new("MetaInfo",
                 data_type = "FLOAT",
                 dims = c(64, 64, 32, 100),
                 spatial_axes = new("AxisSet3D"),
                 additional_axes = new("AxisSet"),
                 spacing = c(3, 3, 4),
                 origin = c(0, 0, 0),
                 label = "fMRI_run1")
```
