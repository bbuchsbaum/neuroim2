# Generic function to position kernel in a position in image space

Generic function to position kernel in a position in image space

## Usage

``` r
embed_kernel(x, sp, center_voxel, ...)

# S4 method for class 'Kernel,NeuroSpace,numeric'
embed_kernel(x, sp, center_voxel, weight = 1)
```

## Arguments

- x:

  the kernel object

- sp:

  the space to embed the kernel

- center_voxel:

  the voxel marking the center of the kernel in the embedded space

- ...:

  extra args

- weight:

  multiply kernel weights by this value

## Value

An object representing the embedded kernel in the specified space.

## Examples

``` r
# Create a 3D Gaussian kernel with dimensions 3x3x3 and voxel size 1x1x1
kern <- Kernel(kerndim = c(3,3,3), vdim = c(1,1,1), FUN = dnorm, sd = 1)

# Create a NeuroSpace object to embed the kernel in
space <- NeuroSpace(c(10,10,10), c(1,1,1))

# Embed the kernel at the center of the space (position 5,5,5)
embedded_kern <- embed_kernel(kern, space, c(5,5,5))

# The result is a SparseNeuroVol with kernel weights centered at (5,5,5)
# We can also scale the kernel weights by using the weight parameter
embedded_kern_scaled <- embed_kernel(kern, space, c(5,5,5), weight = 2)

# The scaled kernel has weights twice as large as the original
max(values(embedded_kern_scaled)) == 2 * max(values(embedded_kern))
#> [1] TRUE
```
