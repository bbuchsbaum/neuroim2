---
title: "Four dimensional neuroimaging data"
date: "2022-10-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Four dimensional neuroimaging data}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---



## Working with neuroimaging time-series data

The `neuroim2` package contains data structures and functions for reading, accessing, and processing 4-dimensional neuroimaging data. 

### Reading a four-dimensional NifTI image with read_vec

Here we read in a 4D image consisting of 5 time points,


```r
      library(purrr)
      library(ggplot2)
      file_name <- system.file("extdata", "global_mask_v5.nii", package="neuroim2")
      vec <- read_vec(file_name)
      dim(vec)
#> [1] 64 64 25  5
```

Now imagine we have a set of 4d images. We can read them in with `read_vec`. (Here we are just using three versions of the same file for the example).



```r
    
      file_name <- system.file("extdata", "global_mask_v5.nii", package="neuroim2")
      vec <- read_vec(c(file_name, file_name, file_name))
      dim(vec)
#> [1] 64 64 25 15
      
      vec2 <- read_vec(rep(file_name, 10))
      vec2
#> NeuroVecSeq 
#>   Type           : NeuroVecSeq 
#>   Dimension      : 64 64 25 50 
#>   Spacing        : 3.5  X  3.5  X  3.7 
#>   Origin         : 112  X  -108  X  -46.2 
#>   Axes           : Right-to-Left Posterior-to-Anterior Inferior-to-Superior 
#>   Coordinate Transform : -3.5 0 0 0 0 3.5 0 0 0 0 3.7 0 112 -108 -46.2 1
```

To extract a subset of volumes we can use the `sub_vector` function:


```r
    
      vec_1_6 <- sub_vector(vec, 1:6)
      dim(vec_1_6)
#> [1] 64 64 25  6
```

### Extracting time-series data using the `series` and `series_roi` functions

To get the time-series at voxel (1,1,1) we can use the `series` function:


```r
      
      series(vec_1_6, 1,1,1)
#> [1] 0 0 0 0 0 0
```

We can extract a 4d region of interest with the `series_roi` as follows:



```r
      file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
      vol <- read_vol(file_name)
      roi <- spherical_roi(vol, c(12,12,12), radius=8)
      rvec1 <- series_roi(vec, roi)
      
      ## or alternatively as a pipeline
      rvec2 <- read_vol(file_name) %>% spherical_roi(c(12,12,12), radius=8) %>% series_roi(vec,.)
      rvec2
#> 
#> 
#> ROIVec 
#>   ncol:            49 
#>   nrow:            15 
#>   Parent Dim:      64 64 25 15 
#>   Voxel Cen. Mass: 12 12 12
      
      ## we can extract the ROI values with the `values` method.
      assertthat::assert_that(all(values(rvec1) == values(rvec2)))
#> [1] TRUE
      assertthat::assert_that(all(coords(rvec1) == coords(rvec2)))
#> [1] TRUE
```

We can also extract an ROI using 1d indices:


```r

r1 <- series_roi(vec, 1:100)
r1
#> 
#> 
#> ROIVec 
#>   ncol:            100 
#>   nrow:            15 
#>   Parent Dim:      64 64 25 15 
#>   Voxel Cen. Mass: 27.46 1.36 1
```

Or we can extract a plain matrix using the `series` function:


```r
r2 <- series(vec, 1:100)
dim(r2)
#> [1]  15 100
```

We can also use coordinate indexing using voxel coordinates. First we load a binary mask with the same spatial dimensions as our NeuroVec:


```r
mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))
```

Now we convert indices to voxels and extract a matrix of values at the specified locations:


```r
vox <- index_to_grid(mask, 1:100)

r3 <- series(vec, vox)
dim(r3)
#> [1]  15 100
```

And the same using `series_roi`:


```r
r4 <- series_roi(vec,vox)
r4
#> 
#> 
#> ROIVec 
#>   ncol:            100 
#>   nrow:            15 
#>   Parent Dim:      64 64 25 15 
#>   Voxel Cen. Mass: 27.46 1.36 1
```







