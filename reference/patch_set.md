# Generate a set of coordinate "patches" of fixed size from an image object.

Generate a set of coordinate "patches" of fixed size from an image
object.

## Usage

``` r
patch_set(x, dims, mask, ...)
```

## Arguments

- x:

  the object to extract patches from

- dims:

  a vector indicating the dimensions of the patches

- mask:

  mask indicating the valid patch area

- ...:

  additional args

## Value

A `list` of coordinate patches, each representing a fixed-size region of
the input object.

## Examples

``` r
# Create a simple 3D volume
space <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
vol <- NeuroVol(array(rnorm(1000), c(10,10,10)), space)

# Create a mask with some active voxels
mask <- LogicalNeuroVol(vol > 0, space)

# Extract 3x3x3 patches centered at each active voxel
patches <- patch_set(vol, dims=c(3,3,3), mask=mask)

# Access the first patch
patch1 <- patches[[1]]
dim(patch1)  # Should be c(27) (flattened 3x3x3 patch)
#> NULL
```
