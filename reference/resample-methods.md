# Resample an Image to Match the Space of Another Image

This function resamples a source image to match the spatial properties
(dimensions, resolution, and orientation) of a target image.

This method resamples a NeuroVol object (`source`) to match the
dimensions and orientation of a NeuroSpace object (`target`).

This method preserves discrete cluster labels and label mappings when
resampling clustered volumes to a new space.

## Usage

``` r
resample(source, target, ...)

# S4 method for class 'NeuroVol,NeuroVol'
resample(source, target, interpolation = 3L)

# S4 method for class 'NeuroVol,NeuroSpace'
resample(source, target, interpolation = 3L)

# S4 method for class 'ClusteredNeuroVol,NeuroSpace'
resample(source, target, interpolation = 0L)

# S4 method for class 'ClusteredNeuroVol,NeuroVol'
resample(source, target, interpolation = 0L)
```

## Arguments

- source:

  A NeuroVol object representing the source volume to be resampled.

- target:

  A NeuroSpace object representing the target space to match the
  dimensions and orientation of the source volume.

- ...:

  Additional arguments passed to the resampling function, such as
  interpolation method, boundary handling, or other resampling options.

- interpolation:

  A single integer specifying the type of interpolation to be applied to
  the final resampled image. May be 0 (nearest neighbor), 1 (trilinear),
  or 3 (cubic spline). No other values are valid.

## Value

An object representing the resampled `source` image, with the same
spatial properties as `target`.

## See also

[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.md)
for the base volume class

## Examples

``` r
img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
rspace <- space(img)


newtrans4X3 <- trans(img)[1:4, 1:3]
newtrans4X3 <- newtrans4X3 * c(.5,.5,.5,1)
newtrans <- cbind(newtrans4X3, c(space(img)@origin,1))

rspace <- NeuroSpace(rspace@dim*2, rspace@spacing/2, origin=rspace@origin, trans=trans(img))
# \donttest{
rvol <- resample(img, rspace)
# }


# Create source and target volumes
src_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
targ_vol <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Resample source to match target
resampled <- resample(src_vol, targ_vol, interpolation=1)

```
