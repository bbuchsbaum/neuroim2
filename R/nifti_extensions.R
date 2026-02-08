#' @include all_class.R
#' @include all_generic.R
NULL

# ============================================================================
# NIfTI Extension Functions
# ============================================================================

#' Create a NIfTI Extension
#'
#' @description
#' Constructor function for creating a \code{\link{NiftiExtension-class}} object
#' with proper padding to ensure the size is a multiple of 16 bytes.
#'
#' @param ecode Integer extension code. See \code{\link{NiftiExtensionCodes}} for
#'   known codes. Common values: 4 (AFNI), 6 (comment), 32 (CIFTI).
#' @param data The extension data. Can be:
#'   \itemize{
#'     \item A character string (will be converted to raw with null terminator)
#'     \item A raw vector (used as-is)
#'   }
#'
#' @return A \code{\link{NiftiExtension-class}} object.
#'
#' @details
#' The function automatically handles padding to ensure the total extension size
#' (esize) is a multiple of 16 bytes, as required by the NIfTI specification.
#' The esize includes the 8-byte header (esize + ecode fields).
#'
#' @examples
#' # Create a comment extension
#' ext <- NiftiExtension(ecode = 6L, data = "This is a comment")
#' ext@ecode
#' ext@esize
#'
#' # Create an AFNI extension with XML data
#' afni_xml <- '<?xml version="1.0"?><AFNI_attributes></AFNI_attributes>'
#' afni_ext <- NiftiExtension(ecode = 4L, data = afni_xml)
#'
#' @seealso
#' \code{\link{NiftiExtension-class}}, \code{\link{NiftiExtensionCodes}}
#'
#' @export
NiftiExtension <- function(ecode, data) {
  # Convert ecode to integer
  ecode <- as.integer(ecode)

 # Convert data to raw vector
  if (is.character(data)) {
    # Add null terminator and convert to raw
    edata <- c(charToRaw(data), as.raw(0))
  } else if (is.raw(data)) {
    edata <- data
  } else {
    stop("'data' must be a character string or raw vector")
  }

  # Calculate required size (must be multiple of 16)
  # esize includes the 8-byte header (esize + ecode)
  min_size <- 8L + length(edata)
  padding_needed <- (16L - (min_size %% 16L)) %% 16L
  esize <- min_size + padding_needed

  # Pad edata with null bytes if needed
  if (padding_needed > 0L) {
    edata <- c(edata, raw(padding_needed))
  }

  new("NiftiExtension",
      ecode = ecode,
      esize = esize,
      edata = edata)
}


#' Read NIfTI Extensions from File
#'
#' @description
#' Reads all extension blocks from a NIfTI file. Extensions are located between
#' byte 352 (after the standard header + extender bytes) and vox_offset.
#'
#' @param conn An open binary connection to the NIfTI file, positioned at byte 348.
#' @param vox_offset The byte offset to image data (from header).
#' @param endian Character string specifying byte order ("little" or "big").
#'
#' @return A \code{\link{NiftiExtensionList-class}} object containing all extensions,
#'   or an empty list if no extensions are present.
#'
#' @details
#' The function first reads the 4-byte extender field at bytes 348-351. If the

#' first byte is non-zero, extensions are present. It then reads each extension
#' block sequentially until vox_offset is reached.
#'
#' @keywords internal
#' @noRd
read_nifti_extensions <- function(conn, vox_offset, endian) {
  # Read the 4-byte extender field (bytes 348-351)
  extender <- readBin(conn, what = "raw", n = 4)

  # Check if extensions are present (first byte != 0)
  if (extender[1] == as.raw(0)) {
    return(new("NiftiExtensionList"))
  }

  # Extensions are present - read them
  extensions <- list()
  current_pos <- 352L

  while (current_pos < vox_offset) {
    # Read esize (4 bytes)
    esize <- readBin(conn, what = "integer", n = 1, size = 4, endian = endian)

    # Validate esize
    if (is.na(esize) || esize < 16 || (esize %% 16) != 0) {
      warning(sprintf("Invalid extension size %s at byte %d. Stopping extension read.",
                     if (is.na(esize)) "NA" else esize, current_pos))
      break
    }

    # Check if extension would exceed vox_offset
    if (current_pos + esize > vox_offset) {
      warning(sprintf("Extension at byte %d (size %d) would exceed vox_offset (%d). Stopping.",
                     current_pos, esize, vox_offset))
      break
    }

    # Read ecode (4 bytes)
    ecode <- readBin(conn, what = "integer", n = 1, size = 4, endian = endian)

    # Read edata (esize - 8 bytes)
    edata_size <- esize - 8L
    edata <- readBin(conn, what = "raw", n = edata_size)

    # Create extension object
    ext <- new("NiftiExtension",
               ecode = as.integer(ecode),
               esize = as.integer(esize),
               edata = edata)

    extensions <- c(extensions, list(ext))
    current_pos <- current_pos + esize
  }

  new("NiftiExtensionList", extensions)
}


#' Write NIfTI Extensions to File
#'
#' @description
#' Writes a list of NIfTI extensions to a file connection.
#'
#' @param conn An open binary connection to write to.
#' @param extensions A \code{\link{NiftiExtensionList-class}} object or NULL.
#' @param endian Character string specifying byte order ("little" or "big").
#'
#' @return Invisibly returns the total number of bytes written (including extender).
#'
#' @details
#' First writes the 4-byte extender field, then each extension block.
#' If extensions is NULL or empty, writes {0,0,0,0} for the extender.
#'
#' @keywords internal
#' @noRd
write_nifti_extensions <- function(conn, extensions, endian) {
  bytes_written <- 0L

  # Write extender field (4 bytes at position 348)
  if (is.null(extensions) || length(extensions) == 0) {
    writeBin(raw(4), conn)
    return(invisible(4L))
  }

  # Extensions present - write extender with first byte = 1
  writeBin(as.raw(c(1L, 0L, 0L, 0L)), conn)
  bytes_written <- 4L

  # Write each extension
  for (ext in extensions) {
    writeBin(ext@esize, conn, size = 4, endian = endian)
    writeBin(ext@ecode, conn, size = 4, endian = endian)
    writeBin(ext@edata, conn)
    bytes_written <- bytes_written + ext@esize
  }

  invisible(bytes_written)
}


#' Calculate Total Extension Size
#'
#' @description
#' Calculates the total byte size of all extensions including the extender field.
#'
#' @param extensions A \code{\link{NiftiExtensionList-class}} object or NULL.
#'
#' @return Integer giving the total size in bytes (minimum 4 for extender field).
#'
#' @keywords internal
#' @noRd
total_extension_size <- function(extensions) {
  if (is.null(extensions) || length(extensions) == 0) {
    return(4L)  # Just the extender field
  }

  # Extender (4 bytes) + sum of all extension sizes
  4L + sum(vapply(extensions, function(ext) ext@esize, FUN.VALUE = integer(1)))
}


#' Get Extension Code Name
#'
#' @description
#' Returns the name associated with a NIfTI extension code.
#'
#' @param ecode Integer extension code.
#'
#' @return Character string with the extension name, or "unknown" if not found.
#'
#' @examples
#' ecode_name(4L)   # Returns "AFNI"
#' ecode_name(6L)   # Returns "comment"
#' ecode_name(999L) # Returns "unknown"
#'
#' @export
ecode_name <- function(ecode) {
  idx <- match(ecode, NiftiExtensionCodes)
  if (is.na(idx)) {
    return("unknown")
  }
  names(NiftiExtensionCodes)[idx]
}


#' Parse NIfTI Extension Data
#'
#' @description
#' Parses the raw data in a NIfTI extension based on its extension code.
#' Provides specialized parsing for known extension types.
#'
#' @param ext A \code{\link{NiftiExtension-class}} object.
#' @param ... Additional arguments passed to type-specific parsers.
#'
#' @return Parsed data in an appropriate format:
#'   \itemize{
#'     \item ecode 4 (AFNI): An XML document (if xml2 available) or character string
#'     \item ecode 6 (comment): Character string
#'     \item Other codes: Raw vector (unchanged)
#'   }
#'
#' @examples
#' # Parse a comment extension
#' ext <- NiftiExtension(ecode = 6L, data = "Test comment")
#' parse_extension(ext)  # Returns "Test comment"
#'
#' @seealso
#' \code{\link{parse_afni_extension}} for AFNI-specific parsing.
#'
#' @export
parse_extension <- function(ext, ...) {
  if (!is(ext, "NiftiExtension")) {
    stop("'ext' must be a NiftiExtension object")
  }

  switch(as.character(ext@ecode),
    "4" = parse_afni_extension(ext, ...),
    "6" = parse_comment_extension(ext),
    {
      # Default: return raw data with a warning for unknown types
      if (!(ext@ecode %in% NiftiExtensionCodes)) {
        warning(sprintf("Unknown extension code %d. Returning raw data.", ext@ecode))
      }
      ext@edata
    }
  )
}


#' Parse Comment Extension
#'
#' @description
#' Parses a comment extension (ecode = 6) to extract the text content.
#'
#' @param ext A \code{\link{NiftiExtension-class}} object with ecode = 6.
#'
#' @return Character string containing the comment text.
#'
#' @keywords internal
#' @noRd
parse_comment_extension <- function(ext) {
  # Remove trailing null bytes and convert to character
  edata <- ext@edata
  # Find first null byte
  null_pos <- which(edata == as.raw(0))
  if (length(null_pos) > 0) {
    edata <- edata[seq_len(null_pos[1] - 1L)]
  }
  rawToChar(edata)
}


#' Parse AFNI Extension
#'
#' @description
#' Parses an AFNI extension (ecode = 4) containing XML-formatted attributes.
#'
#' @param ext A \code{\link{NiftiExtension-class}} object with ecode = 4.
#' @param as_xml Logical; if TRUE (default) and xml2 is available, returns an
#'   xml_document object. Otherwise returns the raw XML string.
#'
#' @return If \code{as_xml = TRUE} and xml2 is available, returns an xml_document.
#'   Otherwise returns a character string containing the XML.
#'
#' @details
#' AFNI stores dataset attributes in an XML format within the NIfTI extension.
#' The XML contains elements like HISTORY_NOTE, volume labels, tagged points,
#' and other AFNI-specific metadata.
#'
#' @examples
#' \dontrun{
#' # Read a NIfTI file with AFNI extension
#' hdr <- read_nifti_header("afni_file.nii")
#' afni_ext <- hdr$extensions[[1]]
#' parsed <- parse_afni_extension(afni_ext)
#' }
#'
#' @seealso
#' \code{\link{get_afni_attribute}} for extracting specific AFNI attributes.
#'
#' @export
parse_afni_extension <- function(ext, as_xml = TRUE) {
  if (ext@ecode != 4L) {
    stop("Extension is not an AFNI extension (ecode != 4)")
  }

  # Get XML string (remove trailing nulls)
  xml_string <- parse_comment_extension(ext)  # Reuse null-stripping logic

  if (as_xml && requireNamespace("xml2", quietly = TRUE)) {
    return(xml2::read_xml(xml_string))
  }

  xml_string
}


#' Get AFNI Attribute from Extension
#'
#' @description
#' Extracts a specific attribute value from a parsed AFNI extension.
#'
#' @param ext A \code{\link{NiftiExtension-class}} object with ecode = 4,
#'   or an xml_document from \code{\link{parse_afni_extension}}.
#' @param attr_name Character string specifying the attribute name to retrieve
#'   (e.g., "HISTORY_NOTE", "BRICK_LABS").
#'
#' @return The attribute value, or NULL if not found. The type depends on the
#'   attribute's ni_type in the XML.
#'
#' @examples
#' \dontrun{
#' # Get the history note from an AFNI extension
#' history <- get_afni_attribute(afni_ext, "HISTORY_NOTE")
#' }
#'
#' @export
get_afni_attribute <- function(ext, attr_name) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required to parse AFNI attributes")
  }

  # Parse if needed
  if (is(ext, "NiftiExtension")) {
    xml_doc <- parse_afni_extension(ext, as_xml = TRUE)
  } else if (inherits(ext, "xml_document") || inherits(ext, "xml_node")) {
    xml_doc <- ext
  } else {
    stop("'ext' must be a NiftiExtension or xml_document")
  }

  # Find the attribute element
  xpath <- sprintf("//AFNI_atr[@atr_name='%s']", attr_name)
  node <- xml2::xml_find_first(xml_doc, xpath)

  if (is.na(node)) {
    return(NULL)
  }

  # Get the value - text content of the node
  value <- xml2::xml_text(node)

  # Try to get type information
  ni_type <- xml2::xml_attr(node, "ni_type")

  # Convert based on type if available
  if (!is.na(ni_type)) {
    if (grepl("float|Float", ni_type, ignore.case = TRUE)) {
      value <- as.numeric(strsplit(trimws(value), "\\s+")[[1]])
    } else if (grepl("int|Int", ni_type, ignore.case = TRUE)) {
      value <- as.integer(strsplit(trimws(value), "\\s+")[[1]])
    }
    # String type returns as-is
  }

  value
}


#' List AFNI Attributes in Extension
#'
#' @description
#' Returns a character vector of all attribute names in an AFNI extension.
#'
#' @param ext A \code{\link{NiftiExtension-class}} object with ecode = 4,
#'   or an xml_document from \code{\link{parse_afni_extension}}.
#'
#' @return Character vector of attribute names.
#'
#' @examples
#' \dontrun{
#' # List all attributes in an AFNI extension
#' attrs <- list_afni_attributes(afni_ext)
#' print(attrs)
#' }
#'
#' @export
list_afni_attributes <- function(ext) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required to parse AFNI attributes")
  }

  # Parse if needed
  if (is(ext, "NiftiExtension")) {
    xml_doc <- parse_afni_extension(ext, as_xml = TRUE)
  } else if (inherits(ext, "xml_document") || inherits(ext, "xml_node")) {
    xml_doc <- ext
  } else {
    stop("'ext' must be a NiftiExtension or xml_document")
  }

  # Find all AFNI_atr elements and get their names
  nodes <- xml2::xml_find_all(xml_doc, "//AFNI_atr")
  xml2::xml_attr(nodes, "atr_name")
}


# ============================================================================
# Show Methods
# ============================================================================

#' @rdname NiftiExtension-class
#' @param object A \code{\linkS4class{NiftiExtension}} object.
#' @export
setMethod("show", "NiftiExtension", function(object) {
  cat("NIfTI Extension\n")
  cat(sprintf("  ecode: %d (%s)\n", object@ecode, ecode_name(object@ecode)))
  cat(sprintf("  esize: %d bytes\n", object@esize))
  cat(sprintf("  data:  %d bytes\n", length(object@edata)))

  # Show preview for known text types
  if (object@ecode %in% c(4L, 6L)) {
    preview <- parse_comment_extension(object)
    if (nchar(preview) > 60) {
      preview <- paste0(substr(preview, 1, 57), "...")
    }
    cat(sprintf("  preview: \"%s\"\n", preview))
  }
  invisible(object)
})


#' @rdname NiftiExtensionList-class
#' @param object A \code{\linkS4class{NiftiExtensionList}} object.
#' @export
setMethod("show", "NiftiExtensionList", function(object) {
  n <- length(object)
  cat(sprintf("NIfTI Extension List (%d extension%s)\n", n, if (n == 1) "" else "s"))

  if (n > 0) {
    total_bytes <- sum(vapply(object, function(x) x@esize, integer(1)))
    cat(sprintf("  Total size: %d bytes\n", total_bytes))

    for (i in seq_along(object)) {
      ext <- object[[i]]
      cat(sprintf("  [%d] ecode=%d (%s), %d bytes\n",
                 i, ext@ecode, ecode_name(ext@ecode), ext@esize))
    }
  }
  invisible(object)
})


# ============================================================================
# Generic Methods for Extensions
# ============================================================================

#' Get Extensions from an Object
#'
#' @description
#' Generic function to retrieve NIfTI extensions from various object types.
#'
#' @param x An object potentially containing extensions.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link{NiftiExtensionList-class}} object, or NULL if no extensions.
#'
#' @export
setGeneric("extensions", function(x, ...) standardGeneric("extensions"))


#' Get Extension by Code
#'
#' @description
#' Retrieve extensions with a specific extension code from a list.
#'
#' @param x A \code{\link{NiftiExtensionList-class}} object.
#' @param ecode Integer extension code to filter by.
#'
#' @return A \code{\link{NiftiExtensionList-class}} containing only extensions
#'   with the specified code.
#'
#' @export
setGeneric("extension", function(x, ecode) standardGeneric("extension"))


#' @rdname extension
#' @export
setMethod("extension", c("NiftiExtensionList", "numeric"), function(x, ecode) {
  ecode <- as.integer(ecode)
  # Convert to list to ensure proper iteration with vapply
  lst <- as.list(x)
  if (length(lst) == 0) {
    return(new("NiftiExtensionList"))
  }
  matches <- vapply(lst, function(ext) ext@ecode == ecode, logical(1))
  new("NiftiExtensionList", lst[matches])
})


#' Check if Extensions are Present
#'
#' @description
#' Tests whether an object has any NIfTI extensions.
#'
#' @param x An object to test.
#'
#' @return Logical indicating whether extensions are present.
#'
#' @export
setGeneric("has_extensions", function(x) standardGeneric("has_extensions"))


#' @rdname has_extensions
#' @export
setMethod("has_extensions", "NiftiExtensionList", function(x) {
  length(x) > 0
})


#' @rdname has_extensions
#' @export
setMethod("has_extensions", "list", function(x) {
  if (is.null(x$extensions)) {
    return(FALSE)
  }
  length(x$extensions) > 0
})
