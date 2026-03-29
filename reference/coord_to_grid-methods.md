# convert n-dimensional real world coordinates to grid coordinates

convert n-dimensional real world coordinates to grid coordinates

## Usage

``` r
coord_to_grid(x, coords)

# S4 method for class 'NeuroSpace,matrix'
coord_to_grid(x, coords)

# S4 method for class 'NeuroSpace,numeric'
coord_to_grid(x, coords)

# S4 method for class 'NeuroVol,matrix'
coord_to_grid(x, coords)

# S4 method for class 'NeuroVol,numeric'
coord_to_grid(x, coords)
```

## Arguments

- x:

  the object

- coords:

  a matrix of real world coordinates

## Value

A numeric `matrix` of grid coordinates.

## Examples

``` r
# Create a simple 3D volume
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
grid <- coord_to_grid(bvol, coords)
world <- grid_to_coord(bvol, grid)
all.equal(coords, world)
#> [1] TRUE
```
