# Generic function to convert voxel coordinates in the reference space (LPI) to native array space.

Generic function to convert voxel coordinates in the reference space
(LPI) to native array space.

## Usage

``` r
grid_to_grid(x, vox)

# S4 method for class 'NeuroSpace,matrix'
grid_to_grid(x, vox)

# S4 method for class 'matrix,matrix'
grid_to_grid(x, vox)
```

## Arguments

- x:

  the object

- vox:

  a matrix of LPI voxel coordinates

## Value

A numeric `matrix` of native voxel coordinates.

## Examples

``` r
# Create a simple 3D volume in LPI orientation
space <- NeuroSpace(c(10,10,10), c(2,2,2))

# Create a reoriented space in RAS orientation
space_ras <- reorient(space, c("R", "A", "S"))

# Convert coordinates between orientations
voxel_coords <- t(matrix(c(1,1,1)))
new_coords <- grid_to_grid(space_ras, voxel_coords)
print(new_coords)
#>      [,1] [,2] [,3]
#> [1,]   10   10   10
```
