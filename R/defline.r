#' @include all-generics.r utils.r blast-utils.r class-utils.r
#' @importFrom stringr str_split_fixed str_extract_all perl
NULL

# Defline-class ----------------------------------------------------------

.valid.Defline <- function(object) {
  errors <- character()
  valid_tag <- object@tag %in% .tags
  if (!all(valid_tag) || is.na(object@tag)) {
    msg <- paste0("Invalid defline tag ", sQuote(object@tag)[!valid_tag])
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0L) TRUE else errors
}

#' Class \code{"Defline"}
#' 
#' @slot tag
#' @slot accession
#' @slot locus
#' @slot description
#' @slot species
#' @keywords internal classes
#' @examples 
#' showClass("Defline")
setClass(Class = "Defline",
         slots = c(tag = "character",
                   accession = "character",
                   locus = "character",
                   description = "character",
                   species = "character"),
         validity = .valid.Defline)


# DeflineSet-class --------------------------------------------------------


#' Class \code{"DeflineSet"} of \code{"Defline"} objects
#' 
#' @slot .Data Inherited from the \code{\link{list}} class.
#' @keywords internal classes
#' @examples 
#' showClass("DeflineSet")
setClass(Class = "DeflineSet", contains = "list",
         validity = listclassValidator('DeflineSet', 'Defline'))

## constructor
DeflineSet <- listclassConstructor('DeflineSet', 'Defline')


# getter, Defline, DeflineSet --------------------------------------------

setMethod(".deflineID", "Defline", function(x, ...) {
  paste0( x@tag %|NA|% '',
          ifelse(is.na(x@tag), '', '|'),
          x@accession %|NA|% '',
          ifelse(is.na(x@locus), '', '|'),
          x@locus %|NA|% '', collapse='|' )
})

setMethod(".deflineID", "DeflineSet", function(x, ...) {
  vapply(x, .deflineID, character(1))
})

setMethod(".deflineDesc", "Defline", function(x, with_species = TRUE) {
  if (with_species) {
    paste0( x@description %|na|% '',
            ifelse(is.na(x@species), ' ', ' ['),
            x@species %|na|% '',
            ifelse(is.na(x@species), '', ']') )
  } else {
    x@description %|na|% ''
  }
})

setMethod(".deflineDesc", "DeflineSet", function(x, with_species = TRUE) {
  vapply(x, .deflineDesc, with_species = with_species, FUN.VALUE = character(1))
})

setMethod(".getDeflineID", 'Defline', function(x, db = 'any') {
  db <- match.arg(db, c('any', .tags))
  if (db == 'any') {
    "Not implemented"
  } else if (db == 'gnl') {
    # General database identifier: gnl|database|identifier
    ans <- x@locus[x@tag == 'gnl']
    attr(ans, 'database') <- x@accession
    ans
  } else {
    x@accession[x@tag == db] %||% NA_character_
  }
})

setMethod(".getDeflineID", 'DeflineSet', function(x, db = 'any') {
  if (db == 'any') {
    "Not implemented"
  } else if (db == 'gnl') {
    ans <- vapply(x, .getDeflineID, db = 'gnl', FUN.VALUE = character(1))
    attr(ans, 'database') <- vapply(x, slot, 'accession', FUN.VALUE = character(1))
    ans
  } else {
    vapply(x, .getDeflineID, db = db, FUN.VALUE = character(1))
  }
})


# show, Defline, DeflineSet ----------------------------------------------


.show_defline <- function(x, show = TRUE) {
  showme <- sprintf("%s %s", .deflineID(x), .deflineDesc(x))
  if (show) {
    cat(showme, sep='\n')
  }
  invisible(showme)
}

setMethod("show", "Defline", function(object) .show_defline(object))

setMethod("show", "DeflineSet", function(object) {
  x <- lapply(object, .show_defline, show = FALSE)
  cat(sprintf('[%s] %s', seq_along(x), x), sep='\n')
  invisible(unlist(x))
})

setAs('Defline', 'character', function(from) {
  def <- .show_defline(from, show = FALSE)
  def
})

setAs('DeflineSet', 'character', function(from) {
  def <- vapply(from, .show_defline, show = FALSE, FUN.VALUE = character(1))
  def
})

# j <- 1
# t <- tags[[j]]
# i <- ids[[j]]
# d <- descriptions[[j]]
# s <- species[[j]]
.parseDeflines <- function(tags, ids, descriptions, species) {  
  DeflineSet(
    Map( function(t, i, d, s) {
      x <- compactChar(unlist(strsplit(i, paste(t, collapse="|"))))
      if (!all(is.na(t))) {
        x <- str_split_fixed(x, '\\|', 3)
        accn <- x[, 2L] %|% NA_character_
        loc <- x[, 3L] %|% NA_character_
      } else {
        accn <- x
        loc <- NA_character_
      }
      new("Defline", tag = t, accession = accn, locus = loc, description = d, species = s)
    }, t = tags, i = ids, d = descriptions, s = species, USE.NAMES = FALSE)
  )
}

#' Parse fasta definition lines
#' 
#' @param x List or character vector of NCBI fasta deflines.
#' @keywords internal
Deflines <- function(x) {
  # first split into identifier and description at the first blank space
  x <- compactChar(unlist(strsplit(unlist(x), "^>")))
  x <- str_split_fixed(x, " ", 2)
  ids <- x[, 1]
  x <- x[, 2]
  x <- str_split_fixed(x, " \\[|\\]", 3)
  descriptions <- as.list(x[, 1] %|% NA_character_)
  species <- as.list(x[, 2] %|% NA_character_)
  
  # parse identifier patterns
  # first we extract the database tags which always are 2 or 3 lowercase
  # letters followed by a pipe.
  db_pattern <- perl("([[:lower:]]{2,3})(?=\\|)")
  tags <- str_extract_all(ids, db_pattern) %|% NA_character_
  ids <- as.list(ids)
  .parseDeflines(tags, ids, descriptions, species)
}
