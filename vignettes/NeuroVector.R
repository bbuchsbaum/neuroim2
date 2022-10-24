## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(purrr)
library(assertthat)
library(neuroim2)
options(mc.cores=1)

## -----------------------------------------------------------------------------
      library(purrr)
      library(ggplot2)
      file_name <- system.file("extdata", "global_mask_v5.nii", package="neuroim2")
      vec <- read_vec(file_name)
      dim(vec)

## -----------------------------------------------------------------------------
    
      file_name <- system.file("extdata", "global_mask_v5.nii", package="neuroim2")
      vec <- read_vec(c(file_name, file_name, file_name))
      dim(vec)
      
      vec2 <- read_vec(rep(file_name, 10))
      vec2

## -----------------------------------------------------------------------------
    
      vec_1_6 <- sub_vector(vec, 1:6)
      dim(vec_1_6)

## -----------------------------------------------------------------------------
      
      series(vec_1_6, 1,1,1)

## -----------------------------------------------------------------------------
      file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
      vol <- read_vol(file_name)
      roi <- spherical_roi(vol, c(12,12,12), radius=8)
      rvec1 <- series_roi(vec, roi)
      
      ## or alternatively as a pipeline
      rvec2 <- read_vol(file_name) %>% spherical_roi(c(12,12,12), radius=8) %>% series_roi(vec,.)
      rvec2
      
      ## we can extract the ROI values with the `values` method.
      assertthat::assert_that(all(values(rvec1) == values(rvec2)))
      assertthat::assert_that(all(coords(rvec1) == coords(rvec2)))
      

## -----------------------------------------------------------------------------

r1 <- series_roi(vec, 1:100)
r1

## -----------------------------------------------------------------------------
r2 <- series(vec, 1:100)
dim(r2)

## -----------------------------------------------------------------------------
mask <- read_vol(system.file("extdata", "global_mask.nii", package="neuroim2"))

## -----------------------------------------------------------------------------
vox <- index_to_grid(mask, 1:100)

r3 <- series(vec, vox)
dim(r3)

## -----------------------------------------------------------------------------
r4 <- series_roi(vec,vox)
r4

