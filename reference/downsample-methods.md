# Downsample an Image

This function downsamples a neuroimaging object, reducing its spatial
resolution while preserving the temporal dimension.

## Usage

``` r
downsample(x, ...)

# S4 method for class 'DenseNeuroVec'
downsample(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box")

# S4 method for class 'SparseNeuroVec'
downsample(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box")

# S4 method for class 'NeuroVec'
downsample(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box")

# S4 method for class 'DenseNeuroVol'
downsample(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box")

# S4 method for class 'NeuroVol'
downsample(x, spacing = NULL, factor = NULL, outdim = NULL, method = "box")
```

## Arguments

- x:

  A DenseNeuroVol object to downsample

- ...:

  Additional arguments passed to specific downsample methods.

- spacing:

  Target voxel spacing (numeric vector of length 3)

- factor:

  Downsampling factor (single value or vector of length 3, between 0 and
  1)

- outdim:

  Target output dimensions (numeric vector of length 3)

- method:

  Downsampling method (currently only "box" for box averaging)

## Value

An object of the same class as `x`, downsampled according to the
specified parameters.

## Examples

``` r
# Create a sample 4D image
data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
space <- NeuroSpace(dim = c(64, 64, 32, 10), 
                    origin = c(0, 0, 0),
                    spacing = c(2, 2, 2))
nvec <- DenseNeuroVec(data, space)

# Downsample by factor
nvec_down1 <- downsample(nvec, factor = 0.5)

# Downsample to target spacing
nvec_down2 <- downsample(nvec, spacing = c(4, 4, 4))

# Downsample to target dimensions
nvec_down3 <- downsample(nvec, outdim = c(32, 32, 16))

# Create a sample 3D volume
data <- array(rnorm(64*64*32), dim = c(64, 64, 32))
space <- NeuroSpace(dim = c(64, 64, 32), 
                    origin = c(0, 0, 0),
                    spacing = c(2, 2, 2))
vol <- DenseNeuroVol(data, space)

# Downsample by factor
vol_down1 <- downsample(vol, factor = 0.5)

# Downsample to target spacing
vol_down2 <- downsample(vol, spacing = c(4, 4, 4))

# Downsample to target dimensions
vol_down3 <- downsample(vol, outdim = c(32, 32, 16))
```
