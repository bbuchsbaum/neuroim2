# Find 3D anatomical orientation from axis abbreviations

Creates a 3D anatomical orientation from axis abbreviations.

## Usage

``` r
findAnatomy3D(axis1 = "L", axis2 = "P", axis3 = "I")
```

## Arguments

- axis1:

  Character string for first axis (default: "L" for Left)

- axis2:

  Character string for second axis (default: "P" for Posterior)

- axis3:

  Character string for third axis (default: "I" for Inferior)

## Value

An AxisSet3D object representing the anatomical orientation

## Examples

``` r
# Create orientation with default LPI axes
orient <- findAnatomy3D()
# Create orientation with custom axes
orient <- findAnatomy3D("R", "A", "S")
```
