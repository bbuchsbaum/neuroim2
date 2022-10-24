## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")

## -----------------------------------------------------------------------------
    library(neuroim2)
    file_name <- system.file("extdata", "global_mask.nii", package="neuroim2")
    vol <- read_vol(file_name)

## -----------------------------------------------------------------------------
    print(vol)

## -----------------------------------------------------------------------------
    class(vol)
    
    is.array(vol)
    
    dim(vol)
    
    vol[1,1,1]
    
    vol[64,64,24]
    

## -----------------------------------------------------------------------------
    
    vol2 <- vol + vol
    sum(vol2) == 2 * sum(vol)
    
    vol3 <- vol2 - 2*vol
    all(vol3 == 0)

## -----------------------------------------------------------------------------
    
    vol2 <- as.logical(vol)
    class(vol2)
    print(vol2[1,1,1])

## -----------------------------------------------------------------------------
    x <- array(0, c(64,64,64))

## -----------------------------------------------------------------------------
    bspace <- NeuroSpace(dim=c(64,64,64), spacing=c(1,1,1))
    vol <- NeuroVol(x, bspace)
    vol

## -----------------------------------------------------------------------------
    vol2 <- NeuroVol((vol+1)*25, space(vol))
    max(vol2)
    space(vol2)
  

## ----eval=FALSE---------------------------------------------------------------
#      write_vol(vol2, "output.nii")
#  
#      ## adding a '.gz' extension results ina gzipped file.
#      write_vol(vol2, "output.nii.gz")

## -----------------------------------------------------------------------------
    temp <- tempfile()
    uri = "http://www.bic.mni.mcgill.ca/~vfonov/nihpd/obj1/nihpd_asym_07.5-13.5_nifti.zip"
    download.file(uri,temp)
    cmd <- paste("unzip ", "-d", dirname(temp), temp)
    system(cmd)
    vol <- read_vol(paste0(dirname(temp), "/nihpd_asym_07.5-13.5_t1w.nii"))
    

## -----------------------------------------------------------------------------
   plot(vol, zlevels=seq(1, dim(vol)[3], by=15))
    

