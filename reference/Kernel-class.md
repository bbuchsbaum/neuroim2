# Kernel

A class representing an image kernel for image processing, such as
convolution or filtering operations in brain images.

## Slots

- `width`:

  A `numeric` value representing the width of the kernel in voxels. The
  width is typically an odd number to maintain symmetry.

- `weights`:

  A `numeric` vector containing the weights associated with each voxel
  in the kernel.

- `voxels`:

  A `matrix` containing the relative voxel coordinates of the kernel.
  Each row represents a voxel coordinate as (x, y, z).

- `coords`:

  A `matrix` containing the relative real-world coordinates of the
  kernel, corresponding to the voxel coordinates.
