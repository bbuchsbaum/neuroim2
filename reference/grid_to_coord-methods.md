# Generic function to convert N-dimensional grid coordinates to real world coordinates

Generic function to convert N-dimensional grid coordinates to real world
coordinates

## Usage

``` r
grid_to_coord(x, coords)

# S4 method for class 'NeuroSpace,matrix'
grid_to_coord(x, coords)

# S4 method for class 'NeuroSpace,numeric'
grid_to_coord(x, coords)

# S4 method for class 'NeuroVol,matrix'
grid_to_coord(x, coords)

# S4 method for class 'NeuroObj,ANY'
grid_to_coord(x, coords)

# S4 method for class 'NeuroSpace,data.frame'
grid_to_coord(x, coords)
```

## Arguments

- x:

  the object: a `NeuroSpace`, or anything that has one, such as a
  `NeuroVol` or `NeuroVec`

- coords:

  grid coordinates: a matrix or data frame with one point per row, or a
  vector giving a single point

## Value

A numeric `matrix` of real-world coordinates, one row per point.

## Examples

``` r
# Create a simple 3D volume
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
grid_coords <- matrix(c(1.5,1.5,1.5, 5.5,5.5,5.5), ncol=3, byrow=TRUE)
world <- grid_to_coord(bvol, grid_coords)
grid <- coord_to_grid(bvol, world)
all.equal(grid_coords, grid)
#> [1] TRUE
```
