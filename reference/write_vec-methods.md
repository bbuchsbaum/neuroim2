# Write a 4d image vector to disk

Write a 4d image vector to disk

## Usage

``` r
write_vec(x, file_name, format, data_type, ...)

# S4 method for class 'NeuroHyperVec,character,missing,missing'
write_vec(x, file_name)

# S4 method for class 'NeuroHyperVec,character,character,missing'
write_vec(x, file_name, format, data_type, ...)

# S4 method for class 'NeuroHyperVec,character,missing,character'
write_vec(x, file_name, data_type)

# S4 method for class 'ROIVec,character,missing,missing'
write_vec(x, file_name)

# S4 method for class 'NeuroVec,character,missing,missing'
write_vec(x, file_name)

# S4 method for class 'NeuroVec,character,character,missing'
write_vec(
  x,
  file_name,
  format,
  nbit = FALSE,
  compression = 5,
  chunk_dim = c(10, 10, 10, dim(x)[4])
)

# S4 method for class 'NeuroVec,character,missing,character'
write_vec(x, file_name, data_type)

# S4 method for class 'ROIVec,character,missing,missing'
write_vec(x, file_name)

# S4 method for class 'NeuroVec,character,missing,missing'
write_vec(x, file_name)

# S4 method for class 'NeuroVec,character,character,missing'
write_vec(
  x,
  file_name,
  format,
  nbit = FALSE,
  compression = 5,
  chunk_dim = c(10, 10, 10, dim(x)[4])
)

# S4 method for class 'NeuroVec,character,missing,character'
write_vec(x, file_name, data_type)
```

## Arguments

- x:

  an image object, typically a `NeuroVec` instance.

- file_name:

  output file name.

- format:

  file format string. Since "NIFTI" is the only currently supported
  format, this parameter can be safely ignored and omitted.

- data_type:

  the numeric data type. If specified should be a `character` vector of:
  "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE". Otherwise output
  format will be inferred from R the datatype of the image.

- ...:

  extra args

- nbit:

  set nbit compression

- compression:

  compression level 1 to 9

- chunk_dim:

  the dimensions of each chunk

## Value

Invisibly returns `NULL` after writing the vector to disk.

## Examples

``` r

bvec <- NeuroVec(array(0, c(10,10,10,10)), NeuroSpace(c(10,10,10,10), c(1,1,1)))
# \donttest{
# Create temporary files
tmp1 <- tempfile(fileext = ".nii")

# Write vectors to temporary files
write_vec(bvec, tmp1)

# Clean up
unlink(tmp1)
# }
```
