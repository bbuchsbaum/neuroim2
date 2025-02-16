#' @keywords internal
#' @importFrom stringr str_split str_trim str_sub
#' @noRd
parseIntAttribute <- function(line){
	as.integer(str_split(str_trim(paste(line, collapse=" ")), "\\s+")[[1]])
}

#' Parse float attribute
#'
#' This function parses one or more lines of text to extract numeric values.
#'
#' @param line A character vector representing one or more lines of text.
#' @return A numeric vector parsed from the input.
#' @keywords internal
#' @noRd
parseFloatAttribute <- function(line){
	as.numeric(str_split(str_trim(paste(line, collapse=" ")), "\\s+")[[1]])
}

#' Parse string attribute
#'
#' This function parses a line (or lines) of text delimited by tildes ('~') to extract string values.
#'
#' @param line A character vector representing one or more lines of text.
#' @return A character vector containing the extracted string values.
#' @keywords internal
#' @noRd
parseStringAttribute <- function(line) {
	
	res <- str_split(line, "~")[[1]]
	if (length(res) > 1) {
		res[1] <- str_sub(res[1], 2)
		res[1:(length(res)-1)]
	} else {
		str_sub(res[1:(length(res)-1)], 2)
	}
}

#' Parse header element
#'
#' This function extracts attribute information from a collection of input lines.
#'
#' @param inputLines A character vector containing lines representing an attribute element.
#' @return A list with elements: type, name, count, and content.
#' @keywords internal
#' @noRd
parseElement <- function(inputLines) {
	atype <- str_trim(str_split(inputLines[[1]], "=")[[1]])[2]
	name <- str_trim(str_split(inputLines[[2]], "=")[[1]])[2]
	count <- str_trim(str_split(inputLines[[3]], "=")[[1]])[2]

	content <- if (atype == "string-attribute") {
		parseStringAttribute(inputLines[4:length(inputLines)])
	} else if (atype == "integer-attribute") {
		parseIntAttribute(inputLines[4:length(inputLines)])
	} else if (atype == "float-attribute") {
		parseFloatAttribute(inputLines[4:length(inputLines)])
	} else {
		stop("unrecognized attribute type ", atype)
	}
	list(type=atype, name=name, count=count, content=content)

}

#' Read an AFNI header
#'
#' This function reads an AFNI header file and returns its parsed contents as a list.
#'
#' @param file_name The name of the AFNI header file (ending in .HEAD).
#' @return A list representation of the AFNI header.
#' @keywords internal
#' @noRd
read_afni_header <- function(file_name) {
	inputLines <- scan(file_name, what=character(), sep="\n", blank.lines.skip=FALSE)
	idx <- which(unlist(lapply(inputLines, function(lin) lin == ""))) + 1
	lastIdx <- length(inputLines)
	attlen <- diff(c(idx, lastIdx+1))
	attpos <- cbind(idx, attlen)

	header <- apply(attpos, 1, function(row) {
		a <- row[[1]]
		b <- row[[2]]
		parseElement(inputLines[a:(a+b-1)])
	})

	names(header) <- unlist(lapply(header, "[[", "name"))
	header

}

#' Parse integer attribute
#'
#' This function parses one or more lines of text to extract integer values.
#'
#' @param line A character vector representing one or more lines of text.
#' @return An integer vector parsed from the input.
#' @keywords internal
#' @importFrom stringr str_split str_trim str_sub
#' @noRd
parseIntAttribute <- function(line){
	as.integer(str_split(str_trim(paste(line, collapse=" ")), "\\s+")[[1]])
}



