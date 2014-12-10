#' @include blastReportDB-class.r
NULL

# blastTable-class -------------------------------------------------------


#' Class \code{"BlastTable"}
#'
#' @description
#' An S4 class that serves as a container for data parsed from
#' NCBI BLAST tabular output.
#' 
#' @slot program <\code{character}>; The BLAST flavour that generated the data.
#' @slot query <\code{character}>; Query definition.
#' @slot reference <\code{character}>; Reference for BLAST.
#' @slot database <\code{character}>; Name of the database.
#' @slot bit_score <\code{numeric}>; The bit score of the hsp.
#' @slot evalue <\code{numeric}>; The expect value.
#' @slot mlog.evalue <\code{numeric}>; 
#' @slot accession <\code{character}>; Accession number.
#' @slot geneid <\code{character}>; NCBI GI number.
#' @slot table <\code{data.frame}>; Hit table.
#' @seealso
#'  The constructor \code{\link{blastTable}}; the BLAST classes
#'  \code{\linkS4class{BlastReport}} and \code{\linkS4class{BlastReportDB}} 
#' @export
setClass(Class = "BlastTable",
         slots = c(program = "character", query = "character",
                   database = "character", bit_score = "numeric",
                   evalue = "numeric", mlog.evalue = "numeric",
                   accession = "character", geneid = "character",
                   table = "data.frame"),
         prototype = prototype(program = NA_character_, query = NA_character_,
                               database = NA_character_, bit_score = NA_real_,
                               evalue = NA_real_, mlog.evalue = NA_real_,
                               accession = NA_character_, geneid = NA_character_,
                               table = data.frame()))

setMethod("show", "BlastTable",
          function (object) {
            cat(sprintf("Query:    %s\nProgram:  %s\nDatabase: %s\n\n",
                        linebreak(object@query, offset=10),
                        sQuote(object@program), sQuote(object@database)))
            print(object@table)
            return(invisible(NULL))
          })


# subsetting-methods, BlastTable ####

#' @describeIn BlastTable Return a column from a BLAST hit table.
setMethod("$", "BlastTable", function(x, name) {
  slot(x, "table")[[name]]
})

#' @describeIn BlastTable Return selected elements from a BLAST hit table as a \code{data.frame}.
setMethod("[", "BlastTable", function (x, i, j, ..., drop = TRUE) {
  slot(x, "table")[i, j,...]
})

#' @describeIn BlastTable Return selected columns from a BLAST hit table.
setMethod("[[", "BlastTable", function(x, i) {
  slot(x, "table")[[i]]
})

#' @describeIn BlastTable Get the column names from a BLAST hit table.
setMethod("names", "BlastTable", function(x) {
  names(slot(x, "table"))
})
