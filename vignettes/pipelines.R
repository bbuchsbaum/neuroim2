## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(purrr)
library(assertthat)
library(neuroim2)
options(mc.cores=1)

## -----------------------------------------------------------------------------
      library(purrr)
      library(ggplot2)
      file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
      vol <- read_vol(file_name)
      mask.idx <- which(vol>0)
      vol[mask.idx] <- runif(length(mask.idx))
      comp <- conn_comp(vol, threshold=.9)
      
      plot(comp$index, zlevels=seq(1,25,by=3), cmap=rainbow(255))

## -----------------------------------------------------------------------------
mvals <- vol %>% split_clusters(comp$index) %>% map_dbl( ~ mean(.))

## -----------------------------------------------------------------------------
sdvol <- vol %>% searchlight(radius=5, eager=TRUE) %>% map_dbl( ~ sd(.)) %>% NeuroVol(space=space(vol), indices=which(vol!=0))
plot(sdvol, cmap=rainbow(255))

## -----------------------------------------------------------------------------
k <- 12
knnfvol <- vol %>% searchlight(radius=12, eager=TRUE) %>% map_dbl(function(x) {
  ind <- order((x[x@center_index] - x)^2)[1:k]
  mean(x[ind])
}) %>% NeuroVol(space=space(vol), indices=which(vol!=0))
plot(knnfvol, cmap=rainbow(255), zlevels=seq(1,25,by=3))

## -----------------------------------------------------------------------------
avgvol <- vol %>% searchlight_coords(radius=12) %>% map_dbl(function(x) {
  vals <- vol[x]
  mean(vals[vals!=0])
}) %>% NeuroVol(space=space(vol), indices=which(vol!=0))
plot(avgvol, cmap=rainbow(255), zlevels=seq(1,25,by=3))

## -----------------------------------------------------------------------------
slice_means <- vol %>% slices %>% map_dbl(~ mean(.))
plot(slice_means, type='l', ylab="mean intensity", xlab="slice number")

## -----------------------------------------------------------------------------
vec <- concat(vol,vol,vol,vol,vol)
mean_vec <- vec %>% vols %>% map_dbl(~ mean(.))
sd_vec <- vec %>% vols %>% map_dbl(~ sd(.))
assert_that(length(mean_vec) == dim(vec)[4])
assert_that(length(sd_vec) == dim(vec)[4])

## -----------------------------------------------------------------------------
vec <- concat(vol,vol,vol,vol,vol)
mean_vol <- vec %>% vectors() %>% map_dbl(~ mean(.)) %>% NeuroVol(., space=space(vol))
assert_that(all(dim(mean_vol) == dim(vol)))

