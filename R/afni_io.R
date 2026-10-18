#' Split one string on a regular expression, keeping empty pieces
#'
#' Base replacement for \code{stringr::str_split(x, pattern)[[1]]}: unlike
#' \code{strsplit()} it keeps a trailing empty piece ("a~b~" gives
#' \code{c("a", "b", "")}) and returns \code{""} for an empty string. Patterns
#' are Perl-compatible.
#' @param x A single string (only the first element is used).
#' @param pattern Regular expression to split on (must not match the empty string).
#' @keywords internal
#' @noRd
.str_split1 <- function(x, pattern) {
  x <- x[1]
  if (is.na(x)) return(NA_character_)
  m <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (m[1] == -1L) return(x)
  starts <- c(1L, m + attr(m, "match.length"))
  ends <- c(m - 1L, nchar(x))
  substring(x, starts, ends)
}

#' Trim leading and trailing (Unicode) whitespace, as stringr::str_trim()
#' @keywords internal
#' @noRd
.str_trim <- function(x) {
  trimws(x, whitespace = "[\\h\\v]")
}

#' @keywords internal
#' @noRd
parseIntAttribute <- function(line){
	as.integer(.str_split1(.str_trim(paste(line, collapse=" ")), "[\\h\\v]+"))
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
	as.numeric(.str_split1(.str_trim(paste(line, collapse=" ")), "[\\h\\v]+"))
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
	
	res <- .str_split1(line[1], "~")
	if (length(res) > 1) {
		res[1] <- substring(res[1], 2)
		res[1:(length(res)-1)]
	} else {
		substring(res[1:(length(res)-1)], 2)
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
	atype <- .str_trim(.str_split1(inputLines[[1]], "="))[2]
	name <- .str_trim(.str_split1(inputLines[[2]], "="))[2]
	count <- .str_trim(.str_split1(inputLines[[3]], "="))[2]

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
#' @noRd
parseIntAttribute <- function(line){
	as.integer(.str_split1(.str_trim(paste(line, collapse=" ")), "[\\h\\v]+"))
}



