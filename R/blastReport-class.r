#' @include Iteration-class.r
NULL

# BlastHeader-class ------------------------------------------------------

#' BlastHeader
#' 
#' A container for blast header information:
#'
#' @slot version Version of BLAST used; \code{"character"}.
#' @slot reference Reference for BLAST; \code{"character"}.
#' @slot database Name of the database; \code{"character"}.
#' @seealso \code{\linkS4class{BlastReport}}, \code{\linkS4class{blastParameters}}
#' @keywords internal
new_BlastHeader <- 
  setClass(Class = "BlastHeader",
           slots = c(version = 'character',
                     reference = 'character',
                     database = 'character'))

setMethod('show', 'BlastHeader',
          function(object) {
            showme <- sprintf("Database Name: %s\tProgram: %s\n",
                              object@database, object@version)
            cat(showme)
          })


# blastParameters-class --------------------------------------------------


#' blastParameters
#'
#' A container for blast parameters and statistics.
#' 
#' @slot program The BLAST flavour that generated the data; \code{"character"}.
#' @slot matrix Name of the matrix (\code{NA} for nucleotide blast); \code{"character"}.
#' @slot expect Cutoff value; \code{"numeric"}.
#' @slot penalties Open and extend penalties; \code{"numeric"}.
#' @slot sc_match Match score for nucleotide-nucleotide comparison; \code{"integer"}.
#' @slot sc_mismatch Mismatch penalty for nucleotide-nucleotide comparison; \code{"integer"}.
#' @slot filter Filter string; \code{"character"}.
#' @slot num_sequences Number of sequences in the database; \code{"character"}.
#' @slot num_letters Number of letters in the database; \code{"character"}.
#' @slot hsp_length Effective HSP length; \code{"numeric"}.
#' @slot effective_space Effective search space; \code{"numeric"}.
#' @slot ka_params kappa, lambda, entropy; \code{"numeric"}.
#' @seealso \code{\linkS4class{BlastReport}}, \code{\linkS4class{BlastHeader}}
#' @keywords internal
new_BlastParameters <- 
  setClass(Class = "BlastParameters",
           slots = c(program = "character", matrix = "character",
                     expect = "numeric", penalties = "numeric",
                     sc_match = "integer", sc_mismatch = "integer",
                     filter = "character", num_sequences = "character",
                     num_letters = "character", hsp_length = "numeric",
                     effective_space = "numeric", ka_params = "numeric"),
           prototype = prototype(penalties = c(open = NA_real_, extend = NA_real_),
                                 ka_params = c(k = NA_real_, lambda = NA_real_, h = NA_real_)))

setMethod('show', 'BlastParameters',
          function(object) {
            fmt.params <- paste0("Search Parameters:\n",
                                 "  Program:                %s\n",
                                 "  Expect value:           %s\n",
                                 "  Substitution matrix:    %s\n",
                                 "  Match/Mismatch scores:  %s,%s\n",
                                 "  Gapcosts (open,extend): %s,%s\n",
                                 "  Filter string:          %s\n")
            params <- sprintf(fmt.params, object@program, object@expect,
                              object@matrix, object@sc_match, object@sc_mismatch,
                              object@penalties['open'], object@penalties['extend'],
                              object@filter)
            fmt.db <- paste0("Database:\n",
                             "  Number of letters:      %s\n",
                             "  Number of sequences:    %s\n")
            db <- sprintf(fmt.db, object@num_letters, object@num_sequences)
            fmt.stat <- paste0("Statistics:\n",
                               "  Lambda:               %s\n",
                               "  K:                    %s\n",
                               "  H:                    %s\n")
            stat <- sprintf(fmt.db, object@ka_params['lambda'],
                            object@ka_params['k'], object@ka_params['h'])
            
            cat(params, db, stat, sep='')
          })


# BlastTable-class ------------------------------------------------------


#' BlastReport
#'
#' @description
#' \code{"BlastReport"} is the top-level container for data parsed from NCBI
#' Blast XML output. It contains the component data classes:
#' \code{\linkS4class{BlastHeader}},
#' \code{\linkS4class{BlastParameters}}, and
#' \code{\linkS4class{IterationList}}.
#' 
#' The \code{\linkS4class{Iteration}} elements store results from individual
#' BLAST queries. Each \code{Iteration} holds a \code{\linkS4class{HitList}}
#' with possibly multiple \code{\linkS4class{Hit}s}, which, in turn, can
#' contain multiple high-scoring pairs (\code{\linkS4class{Hsp}s}).
#' 
#' Iterations, Hits, and Hsps can be extracted using the accessors
#'  \code{\link{getIteration}}, \code{\link{getHit}}, and \code{\link{getHsp}},
#'  or by directly subsetting a \code{blastReport} object.
#'  
#' E.g. \code{report[[1]][[1]]} will return the first hit in the first query.
#'  
#' @slot header Header information; \code{\linkS4class{BlastHeader}}.
#' @slot params Blast parameters and statistics \code{\linkS4class{BlastParameters}}.
#' @slot iterations Iterations; \code{\linkS4class{IterationList}}.
#' @seealso
#'  The constructor \code{\link{blastReport}}; the BLAST classes
#'  \code{\linkS4class{blastReportDB}} and \code{\linkS4class{blastTable}} 
#' @export
new_BlastReport <- 
  setClass(Class = "BlastReport",
           slots = c(header = "BlastHeader",
                     parameters = "BlastParameters",
                     iterations = "IterationList")
  )


# getter, BlastReport ----------------------------------------------------


## @return BlastHeader
setMethod("getHeader", "BlastReport", function(x, ...) x@header)

## @return BlastParameters
setMethod("getParams", "BlastReport", function(x) x@parameters)

## @return Iteration|IterationList
setMethod("getIteration", "BlastReport", function(x, drop = TRUE) {
  it <- x@iterations
  if (drop && length(it) == 1)
    it[[1]]
  else
    it
})

## @return list<HitLists>
setMethod("getHit", "BlastReport",
          function(x, n = NULL, drop = TRUE) {
  getHit(getIteration(x), n = n, drop = drop)
})

## @return vector<integer>
setMethod("getIterNum", "BlastReport", function(x) {
  getIterNum(getIteration(x))
})

## @return vector<character>
setMethod("getQueryID", "BlastReport", function(x) {
  getQueryID(getIteration(x))
})

## @return vector<character>
setMethod("getQueryDef", "BlastReport", function(x) {
  getQueryDef(getIteration(x))
})

## @return vector<integer>
setMethod("getQueryLen", "BlastReport", function(x) {
  getQueryLen(getIteration(x))
})


# subsetting, blastReport ------------------------------------------------


## Subset to IterationList
setMethod("[", "BlastReport", function(x, i, j, ..., drop) {
  x@iterations[i]
})

## Subset to Iteration
setMethod("[[", "BlastReport", function(x, i, j, ...) {
  x@iterations[[i]]
})

setMethod("is.na", "BlastReport", function(x) {
  vapply(x@iterations, is.na, FALSE, USE.NAMES=FALSE)
})


# show, blastReport ------------------------------------------------------


.show_BlastReport <- function(object) {
  olen <- length(object@iterations)
  cat(sprintf("A %s instance with %s iteration%s.\n",
              sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
      sep="")
  show( getHeader(object) )
  cat('\n')
  op <- options("showHits" = getOption("showHits", default = 3L))
  x <- lapply(object@iterations, .show_Iteration)
  options(op)
}


#' @details
#' The \code{show} methods for various blast objects can be modified by
#' a number of global options. Specifically, the number of
#' \code{\linkS4class{Hit}s} shown when displaying
#' \code{\linkS4class{Iteration}s} or \code{\linkS4class{HitList}s} are set
#' by \code{showHits}. Whether alignments are displayed or not is controlled
#' by \code{showAlignment}. Setting these options to \code{NULL}, restores
#' the defaults.
#'  
#' See the general documentation of \code{\link[methods]{show}} method for
#' the expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' @describeIn BlastReport
#' @export
#' @examples
#' # options("showHits" = 20)
#' # show(hitlist)
#' # options("showHits" = NULL)
#' # show(hitlist)
setMethod("show", "BlastReport",
          function(object) {
            .show_BlastReport(object)
          })
