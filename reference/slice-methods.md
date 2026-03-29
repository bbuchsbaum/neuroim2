# Extract image slice

Extract a 2D slice from an image volume

## Usage

``` r
slice(x, zlevel, along, orientation, ...)

# S4 method for class 'NeuroVol,numeric,numeric,missing'
slice(x, zlevel, along, orientation)

# S4 method for class 'NeuroVol,numeric,NeuroSpace,AxisSet3D'
slice(x, zlevel, along, orientation)
```

## Arguments

- x:

  the object

- zlevel:

  coordinate (in voxel units) along the sliced axis

- along:

  the axis along which to slice

- orientation:

  the target orientation of the 2D slice

- ...:

  additional arguments

## Value

A 2D slice from the image volume.
