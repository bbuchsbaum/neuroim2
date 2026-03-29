# Apply a function to an object.

This function applies a function to an object, with additional arguments
passed to the function using the `...` argument. The mapping object
specifies how the function is to be applied, and can take many different
forms, depending on the object and function used. The return value
depends on the function used.

## Usage

``` r
mapf(x, m, ...)

# S4 method for class 'NeuroVol,Kernel'
mapf(x, m, mask = NULL)
```

## Arguments

- x:

  the object that is mapped.

- m:

  the mapping object.

- ...:

  additional arguments to be passed to the function.

- mask:

  restrict application of kernel to masked area

## Value

The result of applying the mapping function to `x`.

## Examples

``` r
# Create a simple 3D volume
bspace <- NeuroSpace(c(10,10,10), c(1,1,1))
vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), bspace)

# Create a 3x3x3 mean smoothing kernel
kern <- Kernel(c(3,3,3),  vdim=c(3,3,3))

# Apply the kernel to smooth the volume
smoothed_vol <- mapf(vol, kern)
```
