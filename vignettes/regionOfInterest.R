## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(assertthat)
library(purrr)
library(neuroim2)
options(mc.cores=1)

## -----------------------------------------------------------------------------
      library(neuroim2)        
      file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
      vol <- read_vol(file_name)

