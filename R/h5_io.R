#library(hdf5r)

# to_nih5_header <- function(vol, h5obj) {
#   hdr <- as_nifti_header(vol, file_name=h5obj$filename)
#   #file.h5 <- H5File$new(test_filename, mode = "w")
#   #file.h5
#
#   hdr_group <- h5obj$create_group("header")
#   for (el in names(hdr)) {
#     hdr_group[[el]] <- hdr[[el]]
#   }
#
#   h5obj
#
# }

to_h5_latentvec <- function(vec, file_name=NULL, data_type="FLOAT",
                            chunk_dim==NULL, nbit=FALSE, compression=6) {

  assert_that(inherits(vec, "LatentNeuroVec"))

  if (!endsWith(file_name, ".lv.h5")) {
    file_name <- paste0(file_name, ".lv.h5")
  }

  h5obj <- hdf5r::H5File$new(file_name)

  hdf5r::h5attr(h5obj, "rtype") <- class(vec)

  dgroup <- h5obj$create_group("data")
  dspace <- h5obj$create_group("space")

  dspace[["dim"]] = dim(vec)
  dspace[["spacing"]] <- spacing(vec)
  dspace[["origin"]] <- origin(space(vec))
  dspace[["trans"]] <- trans(space(vec))

  dtype <- hdf5r::H5P_DATASET_CREATE$new()

  if (nbit && compression > 0) {
    dtype <- dtype$set_nbit()
  }


  h5dtype <- switch(data_type,
                    "BINARY"=hdf5r::h5types$H5T_NATIVE_HBOOL,
                    "SHORT"=hdf5r::h5types$H5T_NATIVE_SHORT,
                    "INT"=hdf5r::h5types$H5T_NATIVE_INT,
                    "INTEGER"=hdf5r::h5types$H5T_NATIVE_INT,
                    "FLOAT"=hdf5r::h5types$H5T_NATIVE_FLOAT,
                    "DOUBLE"=hdf5r::h5types$H5T_NATIVE_DOUBLE,
                    "LONG"=hdf5r::h5types$H5T_NATIVE_LONG,
                    NULL)


  if (!is.null(chunk_dim)) {
    dtype <- dtype$set_chunk(chunk_dim)$set_fill_value(h5dtype, 0)$set_deflate(compression)
  } else {
    dtype <- dtype$set_fill_value(h5dtype, 0)$set_deflate(compression)
  }

  basis_ds <- hdf5r::H5S$new(dims = dim(vec@basis), maxdims = dim(vec@basis))
  loadings_ds <- hdf5r::H5S$new(dims = dim(vec@loadings), maxdims = dim(vec@loadings))

  basis_dset <- dgroup$create_dataset(name = "basis", space = basis_ds,
                                dtype = h5dtype, dataset_create_pl = dtype,
                                gzip_level = compression)


  loadings_dset <- dgroup$create_dataset(name = "loadings", space = loadings_ds,
                                      dtype = h5dtype, dataset_create_pl = dtype,
                                      gzip_level = compression)




  basis_dset[,] <- as.matrix(vec@basis)
  loadings_dset[,] <- as.matrix(vec@loadings)
  dgroup[["offset"]] <- vec@offset
  dgroup[["indices"]] <- as.integer(vec@map@indices)

  h5obj


}


#' @keywords internal
to_nih5_vec <- function(vec, file_name=NULL, data_type="FLOAT", chunk_dim=c(4,4,4,dim(vec)[4]),
                        nbit=FALSE, compression=6) {

  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package \"hdf5r\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  vec <- as(vec, "DenseNeuroVec")

  assert_that(compression >=0 && compression <=9)
  if (is.null(file_name)) {
    file_name <- paste0(tempfile(), ".h5")
  }

  h5obj <- hdf5r::H5File$new(file_name)

  space_ds <- hdf5r::H5S$new(dims = dim(vec), maxdims = dim(vec))
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


  dtype <- dtype$set_chunk(chunk_dim)$set_fill_value(h5dtype, 0)$set_deflate(compression)
  if (nbit && compression > 0) {
    dtype <- dtype$set_nbit()
  }

  hdf5r::h5attr(h5obj, "rtype") <- class(vec)
  #h5obj[["type"]] <- class(vec)

  dgroup <- h5obj$create_group("data")
  dspace <- h5obj$create_group("space")
  dset <- dgroup$create_dataset(name = "elements", space = space_ds,
                                dtype = h5dtype, dataset_create_pl = dtype, chunk_dims = chunk_dim,
                                gzip_level = compression)

  dspace[["dim"]] = dim(vec)
  dspace[["spacing"]] <- spacing(vec)
  dspace[["origin"]] <- origin(space(vec))
  dspace[["trans"]] <- trans(space(vec))

  dset[,,,] <- as.array(vec)
  h5obj

}

#' @keywords internal
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
  #weather_ds_type$to_text()

  dgroup <- h5obj$create_group("data")
  dset <- dgroup$create_dataset(name = "elements", space = space_ds,
                                            dtype = h5dtype,
                                            dataset_create_pl = dtype,
                                            chunk_dims = cdim,
                                            gzip_level = 5)
  dset[,,] <- vol@.Data
  new("H5NeuroVol", h5obj=h5obj)

}
