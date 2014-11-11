#' @include blastReportDB-class.r
NULL

# blastTable-class -------------------------------------------------------


#' BlastTable
#' 
#' BlastTable is an S4 class that provides a container for data retrived
#' by calls to the NCBI Blast utility.
#' 
#' @slot program The BLAST flavour that generated the data; \code{"character"}.
#' @slot query Query definition; \code{"character"}.
#' @slot reference Reference for BLAST; \code{"character"}.
#' @slot database Name of the database; \code{"character"}.
#' @slot bit_score The bit score of the hsp; \code{"numeric"}.
#' @slot evalue The expect value; \code{"numeric"}.
#' @slot mlog.evalue
#' @slot accession Accession number; \code{"character"}.
#' @slot geneid Accession number; \code{"character"}.
#' @slot table Hit table; \code{"data.frame"}.
#' @seealso
#'  The constructor \code{\link{blastTable}}; the BLAST classes
#'  \code{\linkS4class{blastReport}} and \code{\linkS4class{blastReportDB}} 
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


# subsetting-methods, blastTable ####

setMethod("$", "BlastTable",
          function(x, name) {
            slot(x, "table")[[name]]
          })

setMethod("[", "BlastTable",
          function (x, i, j, ..., drop = TRUE) {
            slot(x, "table")[i,j,...]
          })

setMethod("[[", "BlastTable",
          function(x, i) {
            slot(x, "table")[[i]]
          })

setMethod("names", "BlastTable",
          function(x) {
            names(slot(x, "table"))
          })
