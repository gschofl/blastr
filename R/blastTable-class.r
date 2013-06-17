#' @include blastReportDB-class.r
NULL

# blastTable-class -------------------------------------------------------


#' blastTable class
#' 
#' blastTable is an S4 class that provides a container for data retrived
#' by calls to the NCBI Blast utility.
#' 
#' blastReport objects have ten slots:
#' \describe{
#'   \item{program}{}
#'   \item{version}{}
#'   \item{reference}{}
#'   \item{db}{}
#'   \item{bit_score}{}
#'   \item{evalue}{}
#'   \item{mlog.evalue}{}
#'   \item{accession}{}
#'   \item{gi}{}
#'   \item{table}{}
#' }
#' 
#' @param ... Slots for \sQuote{blastTable} instances.
#' @name blastTable-class
#' @rdname blastTable-class
#' @exportClass blastTable
setClass("blastTable",
         representation(program = "character",
                        query = "character",
                        db = "character",
                        bitscore = "numeric",
                        evalue = "numeric",
                        mlog.evalue = "numeric",
                        accession = "character",
                        gi = "character",
                        table = "data.frame"),
         prototype(program = NA_character_,
                   query = NA_character_,
                   db = NA_character_,
                   bitscore = NA_real_,
                   evalue = NA_real_,
                   mlog.evalue = NA_real_,
                   accession = NA_character_,
                   gi = NA_character_,
                   table = data.frame()))

#' @aliases show,blastTable-method
#' @rdname show-methods
setMethod("show", "blastTable",
          function (object) {
            cat(sprintf("Query:    %s\nProgram:  %s\nDatabase: %s\n\n",
                        linebreak(object@query, offset=10),
                        sQuote(object@program), sQuote(object@db)))
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
