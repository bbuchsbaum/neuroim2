The `neuroim2` package is currently in development. The goal of the project is to to provide basic functionality for working with neuroimaging data in R. `neuroim2` It is an effort to modify and upgrade the `neuroim` package (<https://github.com/bbuchsbaum/neuroim>).

Installation
------------

You can install the development version of `neuroim2` (v. `0.1.0`) from Github with:

``` r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/neuroim2")
```

Usage
-----

To read in a volumetric NIFTI formatted image, use the `read_vol` function:

``` r
library(neuroim2)
```

    ## Loading required package: Matrix

``` r
fname <- system.file("extdata", "global_mask.nii", package="neuroim2")
vol <- read_vol(fname)
```

Now, `vol` can be treated as a three-dimensional array:

``` r
v1 <- vol[1,1,1]
vol2 <- vol + vol
vol3 <- vol2 - vol
all(vol == vol3)
```

    ## [1] TRUE

We can create a 4D image, by concatenating several 3D volumes:

``` r
vec <- concat(vol, vol, vol2)
```

    ## [1] "concat4d"

``` r
series1 <- vec[1,1,1,]
length(series1)
```

    ## [1] 3

Vignettes
---------

See examples of use of `neuroim2` in the [vignettes](https://bbuchsbaum.github.io/neuroim2/articles/index.html).
