# Write a 3d image volume to disk

Write a 3d image volume to disk

## Usage

``` r
write_vol(x, file_name, format, data_type, ...)

# S4 method for class 'NeuroVol,character,missing,missing'
write_vol(x, file_name, format, data_type, ...)

# S4 method for class 'ClusteredNeuroVol,character,missing,missing'
write_vol(x, file_name, format, data_type, ...)

# S4 method for class 'NeuroVol,character,character,missing'
write_vol(x, file_name, format, data_type, ...)

# S4 method for class 'ROIVol,character,character,missing'
write_vol(x, file_name, format, data_type, ...)

# S4 method for class 'NeuroVol,character,missing,character'
write_vol(x, file_name, format, data_type, ...)
```

## Arguments

- x:

  an image object, typically a
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  instance.

- file_name:

  output file name

- format:

  file format string. Since "NIFTI" is the only currently supported
  format, this parameter can be safely ignored and omitted.

- data_type:

  output data type, If specified should be one of "UBYTE", "BYTE",
  "SHORT", "USHORT", "INT", "UINT", "LONG", "ULONG", "FLOAT" or
  "DOUBLE". Defaults to the datatype of the file the image was read from
  when its values are still exactly representable in it, and to "FLOAT"
  otherwise; integer targets are fitted with `scl_slope` / `scl_inter`
  rather than truncated. Otherwise output format will be inferred from R
  the datatype of the image.

- ...:

  additional arguments passed to the format-specific writer.

## Value

Invisibly returns `NULL` after writing the volume to disk.

## Details

The output format will be inferred from file extension.

The output format will be inferred from file extension.
`write_vol(x, "out.nii")` outputs a NIFTI file.
`write_vol(x, "out.nii.gz")` outputs a gzipped NIFTI file.

No other file output formats are currently supported.

## Examples

``` r

bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
# \donttest{
tmp1 <- tempfile(fileext = ".nii")
write_vol(bvol, tmp1)
unlink(tmp1)
# }
```
