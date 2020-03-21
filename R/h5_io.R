#library(hdf5r)

to_nih5_header <- function(vol, h5obj) {
  hdr <- as_nifti_header(vol, file_name=h5obj$filename)
  #file.h5 <- H5File$new(test_filename, mode = "w")
  #file.h5

  hdr_group <- h5obj$create_group("header")
  for (el in names(hdr)) {
    hdr_group[[el]] <- hdr[[el]]
  }

  h5obj

}


to_nih5_vol <- function(vol, file_name=NULL, data_type="FLOAT") {
  if (is.null(file_name)) {
    file_name <- paste0(tempfile(), ".h5")
  }

  h5obj <- hdf5r::H5File$new(file_name)

  space_ds <- hdf5r::H5S$new(dims = dim(vol), maxdims = dim(vol))
  dtype <- hdf5r::H5P_DATASET_CREATE$new()

  h5dtype <- switch(data_type,
                    "BINARY"=hdf5r::h5types$H5T_NATIVE_HBOOL,
                    "SHORT"=hdf5r::h5types$H5T_NATIVE_SHORT,
                    "INT"=hdf5r::h5types$H5T_NATIVE_INT,
                    "INTEGER"=hdf5r::h5types$H5T_NATIVE_INT,
                    "FLOAT"=hdf5r::h5types$H5T_NATIVE_FLOAT,
                    "DOUBLE"=hdf5r::h5types$H5T_NATIVE_DOUBLE,
                    "LONG"=hdf5r::h5types$H5T_NATIVE_LONG,
                    NULL)
  if (is.null(h5dtype)) {
    stop(paste("unsupported data_type:", data_type))
  }

  cdim <- c(dim(vol)[1:2],1)

  dtype <- dtype$set_chunk(cdim)$set_fill_value(h5dtype, 0)

  dgroup <- h5obj$create_group("data")
  dset <- dgroup$create_dataset(name = "elements", space = space_ds,
                                            dtype = h5dtype, dataset_create_pl = dtype, chunk_dims = cdim,
                                            gzip_level = 5)
  dset[,,] <- vol@.Data
  new("H5NeuroVol", h5obj=h5obj)

}
