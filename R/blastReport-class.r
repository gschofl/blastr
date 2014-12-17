#' @include Iteration-class.r
NULL

# BlastHeader-class ------------------------------------------------------

#' Class \code{"BlastHeader"}
#' 
#' An S4 class that that serves as a container for BLAST header information:
#'
#' @slot version <\code{character}>; version of BLAST used.
#' @slot reference <\code{character}>; reference for BLAST.
#' @slot database <\code{character}>; name of the database.
#' @seealso \code{\linkS4class{BlastReport}}, \code{\linkS4class{BlastParameters}}
#' @keywords internal
#' @examples 
#' showClass("BlastHeader")
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


# BlastParameters-class --------------------------------------------------


#' Class \code{"BlastParameters"}
#'
#' @description
#' An S4 class that that serves as the container for blast parameters and
#' statistics.
#' 
#' @slot program \code{character} vector; the BLAST flavour that generated the
#'   data.
#' @slot matrix \code{character} vector; name of the matrix (\code{NA} for
#'   nucleotide blast).
#' @slot expect \code{numeric} vector; cutoff value.
#' @slot penalties \code{numeric} vector; open and extend penalties.
#' @slot sc_match \code{integer} vector; match score for nucleotide-nucleotide
#'   comparison.
#' @slot sc_mismatch \code{integer} vector; mismatch penalty for nucleotide-
#'   nucleotide comparison.
#' @slot filter \code{character} vector; filter string.
#' @slot num_sequences \code{character} vector; number of sequences in the 
#'   database.
#' @slot num_letters \code{character} vector; number of letters in the database.
#' @slot hsp_length \code{numeric} vector; effective HSP length.
#' @slot effective_space \code{numeric} vector; effective search space.
#' @slot ka_params \code{numeric} vector; kappa, lambda, entropy.
#' @seealso \code{\linkS4class{BlastReport}}, \code{\linkS4class{BlastHeader}}
#' @keywords internal
#' @examples 
#' showClass("BlastParameters")
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
                               "  Lambda:                 %s\n",
                               "  K:                      %s\n",
                               "  H:                      %s\n")
            stat <- sprintf(fmt.stat, object@ka_params['lambda'],
                            object@ka_params['k'], object@ka_params['h'])
            
            cat(params, db, stat, sep = '')
          })


# BlastReport-class ------------------------------------------------------


#' Class \code{"BlastReport"}
#'
#' @description
#' An S4 class that that serves as the top-level container for data parsed from
#' NCBI BLAST XML output. It contains the following top-level components:
#' 
#' \itemize{
#'    \item \code{\linkS4class{BlastHeader}}
#'    \item \code{\linkS4class{BlastParameters}}
#'    \item \code{\linkS4class{IterationList}}
#' }
#' 
#' @details
#' The \code{\linkS4class{Iteration}} elements store results from individual
#' BLAST queries. Each \code{Iteration} holds a \code{\linkS4class{HitList}}
#' with possibly multiple \code{\linkS4class{Hit}s}, which, in turn, can
#' contain multiple high-scoring pairs (\code{\linkS4class{Hsp}s}).
#' 
#' Iterations, Hits, and Hsps can be extracted using the accessors
#' \code{\link{getIteration}}, \code{\link{getHit}}, and \code{\link{getHsp}},
#'  or by directly subsetting a \code{BlastReport} object.
#'  
#' E.g. \code{report[[1]][[1]]} will return the first hit in the first query.
#'  
#' @slot header Header information; \code{\linkS4class{BlastHeader}}.
#' @slot params Blast parameters and statistics \code{\linkS4class{BlastParameters}}.
#' @slot iterations Iterations; \code{\linkS4class{IterationList}}.
#' @seealso
#'  The constructor \code{\link{blastReport}}; the BLAST classes
#'  \code{\linkS4class{BlastReportDB}} and \code{\linkS4class{BlastTable}}
#' @keywords classes 
#' @export
new_BlastReport <- 
  setClass(Class = "BlastReport",
           slots = c(header = "BlastHeader",
                     parameters = "BlastParameters",
                     iterations = "IterationList")
  )


# getter, BlastReport ----------------------------------------------------


#' @describeIn BlastReport Return the \code{\linkS4class{BlastHeader}}.
setMethod("getHeader", "BlastReport", function(x, ...) x@header)

#' @describeIn BlastReport Return the \code{\linkS4class{BlastParameters}}.
setMethod("getParams", "BlastReport", function(x) x@parameters)

#' @describeIn BlastReport Return the \code{\linkS4class{Iteration}} or 
#'   \code{\linkS4class{IterationList}}.
setMethod("getIteration", "BlastReport", function(x, i, drop = TRUE) {
  it <- if (missing(i)) x@iterations[] else x@iterations[i]
  if (drop && length(it) == 1) {
    it[[1]]
  } else it
})

#' @describeIn BlastReport Return a list of \code{\linkS4class{HitList}}s.
setMethod("getHit", "BlastReport", function(x, i, drop = TRUE) {
  f <- if (missing(i)) getHit else Partial(getHit, i = i)
  lapply(getIteration(x), f, drop = drop)
})

#' @describeIn BlastReport Return iteration numbers; <\code{integer}>.
setMethod("getIterNum", "BlastReport", function(x) {
  getIterNum(getIteration(x))
})

#' @describeIn BlastReport Return query IDs; <\code{character}>.
setMethod("getQueryID", "BlastReport", function(x) {
  getQueryID(getIteration(x))
})

#' @describeIn BlastReport Return query definitions; <\code{character}>.
setMethod("getQueryDef", "BlastReport", function(x) {
  getQueryDef(getIteration(x))
})

#' @describeIn BlastReport Return query lengths; <\code{integer}>.
setMethod("getQueryLen", "BlastReport", function(x) {
  getQueryLen(getIteration(x))
})


# subsetting, blastReport ------------------------------------------------


#' @describeIn BlastReport Subset to return an \code{\linkS4class{IterationList}}.
setMethod("[", "BlastReport", function(x, i, j, ..., drop) {
  if (missing(i)) x@iterations[] else x@iterations[i]
})

#' @describeIn BlastReport Subset to return an \code{\linkS4class{Iteration}}.
setMethod("[[", "BlastReport", function(x, i, j, ...) {
  x@iterations[[i]]
})

#' @describeIn BlastReport Indicate missing \code{\linkS4class{Iteration}}s.
setMethod("is.na", "BlastReport", function(x) {
  vapply(x@iterations, is.na, FALSE, USE.NAMES = FALSE)
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
setMethod("show", "BlastReport", function(object) {
  .show_BlastReport(object)
})
