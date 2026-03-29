# Resample an image with readable method names

A convenience front-end to \[resample()\] that accepts human-friendly
method names and an engine switch. Internally delegates to the S4
\`resample(source, target, interpolation = 0/1/3)\` methods.

## Usage

``` r
resample_to(
  source,
  target,
  method = c("nearest", "linear", "cubic"),
  engine = c("internal"),
  ...
)
```

## Arguments

- source:

  A \`NeuroVol\` (source image)

- target:

  A \`NeuroVol\` or \`NeuroSpace\` to match

- method:

  Interpolation method: \`"nearest"\`, \`"linear"\`, or \`"cubic"\`

- engine:

  Resampling engine. For now only \`"internal"\` is supported.

- ...:

  Reserved for future options

## Value

A \`NeuroVol\` in the target space

## Examples

``` r
# \donttest{
img <- read_vol(system.file("extdata","global_mask_v4.nii", package="neuroim2"))
sp  <- space(img)
sp2 <- NeuroSpace(sp@dim*2, sp@spacing/2, origin=sp@origin, trans=trans(img))
r1  <- resample_to(img, sp2, method = "linear")
# }
```
