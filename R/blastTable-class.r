#' @include blastReportDB-class.r
NULL

# blastTable-class -------------------------------------------------------


#' blastTable-class
#' 
#' blastTable is an S4 class that provides a container for data retrived
#' by calls to the NCBI Blast utility.
#' 
#' blastReport objects have ten slots:
#' \describe{
#'   \item{\code{program}:}{The BLAST flavour that generated the data; \code{"character"}.}
#'   \item{\code{query}:}{Query definition; \code{"character"}.}
#'   \item{\code{reference}:}{Reference for BLAST; \code{"character"}.}
#'   \item{\code{database}:}{Name of the database; \code{"character"}.}
#'   \item{\code{bit_score}:}{The bit score of the hsp; \code{"numeric"}.}
#'   \item{\code{evalue}:}{The expect value; \code{"numeric"}.}
#'   \item{\code{mlog.evalue}:}{}
#'   \item{\code{accession}:}{Accession number; \code{"character"}.}
#'   \item{\code{geneid}:}{Accession number; \code{"character"}.}
#'   \item{\code{table}:}{Hit table; \code{"data.frame"}.}
#' }
#'
#' @seealso
#'  The constructor \code{\link{blastTable}}; the BLAST classes
#'  \code{\linkS4class{blastReport}} and \code{\linkS4class{blastReportDB}} 
#' @name blastTable-class
#' @rdname blastTable-class
#' @exportClass blastTable
setClass("blastTable",
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

#' @aliases show,blastTable-method
#' @rdname show-methods
setMethod("show", "blastTable",
          function (object) {
            cat(sprintf("Query:    %s\nProgram:  %s\nDatabase: %s\n\n",
                        linebreak(object@query, offset=10),
                        sQuote(object@program), sQuote(object@database)))
            print(object@table)
            return(invisible(NULL))
          })


# subsetting-methods, blastTable ####


setMethod("$", "blastTable",
          function(x, name) {
            slot(x, "table")[[name]]
          })


setMethod("[", "blastTable",
          function (x, i, j, ..., drop = TRUE) {
            slot(x, "table")[i,j,...]
          })


setMethod("[[", "blastTable",
          function(x, i) {
            slot(x, "table")[[i]]
          })


setMethod("names", "blastTable",
          function(x) {
            names(slot(x, "table"))
          })
