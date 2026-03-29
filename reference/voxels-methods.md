# extract voxel coordinates

extract voxel coordinates

## Usage

``` r
voxels(x, ...)

# S4 method for class 'Kernel'
voxels(x, center_voxel = NULL)
```

## Arguments

- x:

  the object to extract voxels from

- ...:

  additional arguments to function

- center_voxel:

  the absolute location of the center of the voxel, default is (0,0,0)

## Value

A `matrix` or `vector` representing voxel coordinates from `x`.

## Examples

``` r
# Create a 3D kernel with dimensions 3x3x3 and voxel size 1x1x1
kern <- Kernel(kerndim = c(3,3,3), vdim = c(1,1,1))

# Get voxel coordinates centered at origin (0,0,0)
vox <- voxels(kern)
# Returns a matrix where each row is a voxel coordinate
# relative to the kernel center

# Get voxel coordinates centered at specific point (5,5,5)
vox_centered <- voxels(kern, center_voxel = c(5,5,5))
# Returns coordinates shifted to be centered at (5,5,5)
```
