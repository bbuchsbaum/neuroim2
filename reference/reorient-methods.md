# Remap the grid-to-world coordinates mapping of an image.

Remap the grid-to-world coordinates mapping of an image.

## Usage

``` r
reorient(x, orient)

# S4 method for class 'NeuroSpace,character'
reorient(x, orient)
```

## Arguments

- x:

  the object

- orient:

  the orientation code indicating the "remapped" axes.

## Value

A reoriented version of `x`.

## Details

When `x` is a `NeuroSpace` object, the `orient` argument should be a
character vector of length 3 specifying the desired anatomical
orientation using single-letter codes. Each letter represents an
anatomical direction:

- First position: "R" (Right) or "L" (Left)

- Second position: "A" (Anterior) or "P" (Posterior)

- Third position: "S" (Superior) or "I" (Inferior)

For example, `c("R", "A", "S")` specifies Right-Anterior-Superior
orientation, while `c("L", "P", "I")` specifies Left-Posterior-Inferior
orientation. The orientation codes determine how the voxel grid
coordinates map to real-world anatomical space.

## Examples

``` r
# Create a NeuroSpace object in LPI (Left-Posterior-Inferior) orientation
space <- NeuroSpace(c(64, 64, 40), c(2, 2, 2))

# Reorient to RAS (Right-Anterior-Superior) orientation
# Use individual axis codes: "R" for Right, "A" for Anterior, "S" for Superior
space_ras <- reorient(space, c("R", "A", "S"))

# The transformation matrix will be updated to reflect the new orientation
# Original and reoriented spaces will have different coordinate mappings
coords <- c(32, 32, 20)
orig_world <- grid_to_coord(space, coords)
new_world <- grid_to_coord(space_ras, coords)
```
