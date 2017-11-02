## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(neuroim2)

## ------------------------------------------------------------------------
    file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
    vol <- read_vol(file_name)

## ------------------------------------------------------------------------
    print(vol)

## ------------------------------------------------------------------------
    class(vol)

    is.array(vol)

    dim(vol)

    vol[1,1,1]

    vol[64,64,24]


## ------------------------------------------------------------------------

    vol2 <- vol + vol
    sum(vol2) == 2 * sum(vol)

    vol3 <- vol2 - 2*vol
    all(vol3 == 0)

## ------------------------------------------------------------------------

    vol2 <- as.logical(vol)
    print(vol2[1,1,1])

## ------------------------------------------------------------------------
    # create an 64X64X64 array of zeros
    x <- array(0, c(64,64,64))

    # create a 'BrainSpace' instance that describes the geometry of the image including, at minimu its dimensions and voxel spacing
    bspace <- NeuroSpace(Dim=c(64,64,64), spacing=c(1,1,1))
    vol <- NeuroVolume(x, bspace)
    vol

## ------------------------------------------------------------------------
    vol2 <- NeuroVol((vol+1)*25, space(vol))
    max(vol2)

    space(vol2)

## ------------------------------------------------------------------------
    write_vol(vol2, "output.nii")

    ## adding a '.gz' extension results ina gzipped file.
    write_vol(vol2, "output.nii.gz")

