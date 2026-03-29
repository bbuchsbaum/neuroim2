# Create a Kernel object from a function of distance from kernel center

This function creates a Kernel object using a kernel function (`FUN`)
that takes the distance from the center of the kernel as its first
argument.

## Usage

``` r
Kernel(kerndim, vdim, FUN = dnorm, ...)
```

## Arguments

- kerndim:

  A numeric vector representing the dimensions in voxels of the kernel.

- vdim:

  A numeric vector representing the dimensions of the voxels in real
  units.

- FUN:

  The kernel function taking its first argument representing the
  distance from the center of the kernel (default: `dnorm`).

- ...:

  Additional parameters to the kernel function, `FUN`.

## Value

A Kernel object with the specified dimensions, voxel dimensions, and
kernel function.

## Examples

``` r
kdim <- c(3, 3, 3)
vdim <- c(1, 1, 1)
k <- Kernel(kerndim = kdim, vdim = vdim, FUN = dnorm, sd = 1)
```
