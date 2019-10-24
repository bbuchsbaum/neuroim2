#' @include all_class.R
NULL
#' @include all_generic.R


#' @rdname file_matches-methods
setMethod(f="file_matches", signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			if (header_file_matches(x,file_name)) {
				file.exists(paste(strip_extension(x, file_name), ".", x@data_extension, sep=""))
			} else if (data_file_matches(x,file_name)) {
				file.exists(paste(strip_extension(x, file_name), ".", x@header_extension, sep=""))
			} else {
				FALSE
			}
		})

#' @rdname header_file_matches-methods
setMethod(f="header_file_matches", signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			regexpr(paste(".*", x@header_extension, "$", sep=""), file_name) > 0
		})

#' @rdname data_file_matches-methods
setMethod(f="data_file_matches", signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			regexpr(paste(".*", x@data_extension, "$", sep=""), file_name) > 0
		})

#' @rdname header_file-methods
setMethod(f="header_file",signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			if (header_file_matches(x, file_name)) {
				file_name
			} else if (data_file_matches(x, file_name)) {
				paste(strip_extension(x, file_name), x@header_extension, sep=".")
			} else {
				stop(paste("could not derive header file name from: ", file_name))
			}
		})

#' @rdname data_file-methods
setMethod(f="data_file",signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			if (data_file_matches(x, file_name)) {
				file_name
			} else if (header_file_matches(x, file_name)) {
				paste(strip_extension(x, file_name), x@data_extension, sep=".")
			} else {
				stop(paste("could not derive data file name from: ", file_name))
			}
		})

#' @rdname strip_extension-methods
setMethod(f="strip_extension",signature=signature(x= "FileFormat", file_name="character"),
		def=function(x, file_name) {
			if (header_file_matches(x, file_name)) {
				ret <- strsplit(file_name, paste(x@header_extension, "$", sep=""))[[1]][1]
				substr(ret, 1, nchar(ret)-1)
			} else if (data_file_matches(x, file_name)) {
				ret <- strsplit(file_name, paste(x@data_extension, "$", sep=""))[[1]][1]
				substr(ret, 1, nchar(ret)-1)
			} else {
				stop("file does not match descriptor: " + x)
			}
		})


.read_meta_info <- function(desc, file_name, read_func, constructor) {

	hfile <- header_file(desc, file_name)
	header <- read_func(hfile)
	header$file_name <- hfile
	constructor(desc, header)
}

#' @rdname read_meta_info-methods
#' @export
setMethod(f="read_meta_info",signature=signature(x= "NIFTIFormat"),
		def=function(x, file_name) {
			.read_meta_info(x, file_name, read_nifti_header, NIFTIMetaInfo)
		})

#' @rdname read_meta_info-methods
#' @export
setMethod(f="read_meta_info",signature=signature(x= "AFNIFormat"),
		def=function(x, file_name) {
			.read_meta_info(x, file_name, read_afni_header, AFNIMetaInfo)

		})


find_descriptor <- function(file_name) {
	if (file_matches(NIFTI, file_name)) NIFTI
	else if (file_matches(NIFTI_GZ, file_name)) NIFTI_GZ
	else if (file_matches(NIFTI_PAIR, file_name)) NIFTI_PAIR
	else if (file_matches(NIFTI_PAIR_GZ, file_name)) NIFTI_PAIR_GZ
	else if (file_matches(AFNI, file_name)) AFNI
	else if (file_matches(AFNI_GZ, file_name)) AFNI_GZ
	else NULL
}

AFNI <- new("AFNIFormat",
		file_format="AFNI",
		header_encoding="raw",
		header_extension="HEAD",
		data_encoding="raw",
		data_extension="BRIK")

AFNI_GZ <- new("AFNIFormat",
		file_format="AFNI",
		header_encoding="gzip",
		header_extension="HEAD",
		data_encoding="gzip",
		data_extension="BRIK.gz")

NIFTI <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="raw",
		header_extension="nii",
		data_encoding="raw",
		data_extension="nii")

NIFTI_GZ <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="gzip",
		header_extension="nii.gz",
		data_encoding="gzip",
		data_extension="nii.gz")

NIFTI_PAIR <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="raw",
		header_extension="hdr",
		data_encoding="raw",
		data_extension="img")

NIFTI_PAIR_GZ <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="gzip",
		header_extension="hdr.gz",
		data_encoding="gzip",
		data_extension="img.gz")


