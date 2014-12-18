#' @include all-generics.r defline.r Hsp-class.r Hit-class.r
NULL

# Iteration-class, IterationList-class -----------------------------------

setClassUnion("HitListOrChar", members = c("HitList", "character"))

#' Class \code{"Iteration"}
#' 
#' An S4 class that holds data parsed from an NCBI BLAST XML \sQuote{Iteration} 
#' element. Each Iteration contains data about a query.
#'
#' @slot iter_num The number of the iteration; \code{"integer"}.
#' @slot query_id Query identifier; \code{"character"}.
#' @slot query_def Query definition; \code{"character"}.
#' @slot query_len Query length; \code{"integer"}.
#' @slot hits List of hits; \code{"\linkS4class{HitList}"}.
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @keywords classes
#' @export
#' @examples 
#' showClass("Iteration")
new_Iteration <- 
  setClass(Class = "Iteration",
           slots = c(iter_num = "integer",
                     query_id = "character",
                     query_def = "character",
                     query_len = "integer",
                     hits = "HitListOrChar",
                     query_env = "environment"),
           prototype = prototype(query_env = new.env(parent = emptyenv()))
  )

#' Class \code{"IterationList"} of \code{"Iteration"} objects
#' 
#' An S4 class that holds a list of \code{\linkS4class{Iteration}} objects.
#'
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @slot .Data Inherited from the \code{\link{list}} class.
#' @seealso \code{\linkS4class{BlastReport}}
#' @keywords classes
#' @export
#' @examples 
#' showClass("IterationList")
setClass(Class    = "IterationList",
         contains = "list",
         validity = listclassValidator('IterationList', 'Iteration')
)

## constructor
IterationList <- listclassConstructor('IterationList', 'Iteration')

# Hsp, nhsps, HspNum ####

## @return list<Hsp|HspList>
setMethod("getHsp", "Iteration", function(x, i, drop = TRUE) {
  lapply(x@hits, getHsp, i = i, drop = drop)
})

## @return vector<numeric>
setMethod('nhsps', 'Iteration', function(x) {
  vapply(x@hits, nhsps, numeric(1))
})

## @return list<vector<numeric>>
setMethod('nhsps', 'IterationList', function(x) {
  lapply(x, nhsps)
})


## Hit, nhits, IterNum, QueryID, QueryDef, QueryLen ####

#' @describeIn Iteration Return a \code{\linkS4class{Hit}} or \code{\linkS4class{HitList}}
setMethod("getHit", "Iteration", function(x, i, drop = TRUE, ...) {
  hit <- if (missing(i)) x@hits[] else x@hits[i]
  if (drop && length(hit) == 1) {
    hit[[1]]
  } else hit
})

#' @describeIn IterationList Returns a list of \code{\linkS4class{Hit}}s or \code{\linkS4class{HitList}}s
setMethod("getHit", "IterationList", function(x, i, drop = TRUE, ...) {
  f <- if (missing(i)) getHit else Partial(getHit, i = i)
  lapply(x, f, drop = drop)
})

#' @describeIn Iteration Returns the number of hits; <\code{numeric}>. 
setMethod('nhits', 'Iteration', function(x) length(x@hits))

#' @describeIn IterationList Returns the numbers of hits; <\code{numeric}>.
setMethod('nhits', 'IterationList', function(x) {
  vapply(x, nhits, 0L)
})

## @return integer
setMethod("getIterNum", "Iteration", function(x) x@iter_num)

## @return vector<integer> 
setMethod("getIterNum", "IterationList", function(x) {
  vapply(x, getIterNum, 0L)
})

## @return character
setMethod("getQueryID", "Iteration", function(x) x@query_id)

## @return vector<character> 
setMethod("getQueryID", "IterationList", function(x) {
  vapply(x, getQueryID, "")
})

## @return character
setMethod("getQueryDef", "Iteration", function(x) x@query_def)

## @return vector<character>
setMethod("getQueryDef", "IterationList", function(x) {
  vapply(x, getQueryDef, "")
})

## @return integer
setMethod("getQueryLen", "Iteration", function(x) x@query_len)

## @return vector<integer>
setMethod("getQueryLen", "IterationList", function(x) {
  vapply(x, getQueryLen, 0L)
})

## HitNum, HitLen, Accession, GeneID ####

## @return vector<integer>
setMethod("getHitNum", "Iteration", function(x) {
  vapply(x@hits, getHitNum, 0L)
})

## @return list<vector<integer>>
setMethod("getHitNum", "IterationList", function(x) {
  lapply(x, getHitNum)
})

## @return vector<integer>
setMethod("getHitLen", "Iteration", function(x) {
  vapply(x@hits, getHitLen, 0L)
})

## @return list<vector<integer>> 
setMethod("getHitLen", "IterationList", function(x) {
  lapply(x, getHitLen)
})

## @return vector<character>
setMethod("getAccession", "Iteration", function(x) {
  vapply(x@hits, getAccession, "")
})

## @return list<vector<character>>
setMethod("getAccession", "IterationList", function(x) {
  lapply(x, getAccession)
})

## @return vector<character>
setMethod("getGeneID", "Iteration", function(x) {
  vapply(x@hits, getGeneID, "")
})

## @return list<vector<character>>
setMethod("getGeneID", "IterationList", function(x) {
  lapply(x, getGeneID)
})

## HitID, HitDef, Defline, PrimaryHitDef, PrimaryDefline ####

## @return list<matrix<character>> or list<vector<integer>>
setMethod("getHitID", "Iteration", function(x, db = 'any') {
  lapply(x@hits, getHitID, db = db)
})

## @return list<vector<character>>
setMethod("getHitDef", "Iteration", function(x) {
  lapply(x@hits, getHitDef)
})

## @return list<vector<character>>
setMethod("getDefline", "Iteration", function(x) {
  lapply(x@hits, getDefline)
})

## @return vector<character>
setMethod("getPrimaryHitDef", "Iteration", function(x) {
  vapply(x@hits, getPrimaryHitDef, "")
})
## @return list<vector<character>>
setMethod("getPrimaryHitDef", "IterationList", function(x) {
  lapply(x, getPrimaryHitDef)
})

## @return vector<character>
setMethod("getPrimaryDefline", "Iteration", function(x) {
  vapply(x@hits, getPrimaryDefline, "")
})

## @return list<vector<character>>
setMethod("getPrimaryDefline", "IterationList", function(x) {
  lapply(x, getPrimaryDefline)
})

## Bitscore, Score, Evalue ####

## list<vector<numeric>>
setMethod("getBitscore", "Iteration",
          function(x, max = FALSE, sum = FALSE) {
            getBitscore(x@hits, max = max, sum = sum)
          })

## @return vector<numeric>
setMethod('getMaxBitscore', 'Iteration', function(x) {
  getMaxBitscore(x@hits)
})

## @return list<vector<numeric>>
setMethod('getMaxBitscore', 'IterationList', function(x) {
  lapply(x, getMaxBitscore)
})

## @return vector<numeric>
setMethod('getTotalBitscore', 'Iteration', function(x) {
  getTotalBitscore(x@hits)
})

## @return list<vector<numeric>>
setMethod('getTotalBitscore', 'IterationList', function(x) {
  lapply(x, getTotalBitscore)
})

## @return list<vector<numeric>>
setMethod('getScore', 'Iteration', function(x, max = FALSE) {
  getScore(x@hits, max = max)
})

## @return list<vector<numeric>> 
setMethod('getEvalue', 'Iteration', function(x, max = FALSE) {
  getEvalue(x@hits, max = max)
})

## Identity, Positive, Gaps, AlignLen ####

## @return list<vector<integer>> 
setMethod('getIdentity', 'Iteration', function(x, max = FALSE) {
  getIdentity(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getPositive', 'Iteration', function(x, max = FALSE) {
  getPositive(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getGaps', 'Iteration', function(x, max = FALSE) {
  getGaps(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getAlignLen', 'Iteration', function(x, max = FALSE) {
  getAlignLen(x@hits, max = max)
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

## @return list<vector<integer>> 
setMethod('getQueryFrom', 'Iteration', function(x, max = FALSE) {
  getQueryFrom(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getQueryTo', 'Iteration', function(x, max = FALSE) {
  getQueryTo(x@hits, max = max)
})

## @return list<vector<integer> 
setMethod('getHitFrom', 'Iteration', function(x, max = FALSE) {
  getHitFrom(x@hits, max = max)
})

## @return list<vector<integer> 
setMethod('getHitTo', 'Iteration', function(x, max = FALSE) {
  getHitTo(x@hits, max = max)
})

## QueryFrame, HitFrame ####

## @return list<vector<integer>>
setMethod('getQueryFrame', 'Iteration', function(x, max = FALSE) {
  getQueryFrame(x@hits, max = max)
})

## @return list<vector<integer>>
setMethod('getHitFrame', 'Iteration', function(x, max = FALSE) {
  getHitFrame(x@hits, max = max)
})

## QueryRange, HitRange ####

## @return list<IRanges> 
setMethod("getQueryRange", "Iteration", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getQueryRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
})

## @return list<IRanges> 
setMethod("getHitRange", "Iteration", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getHitRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
}) 

## getQuerySeq, getHitSeq, getMatch ####

## @return list<BStringSet>
setMethod('getQuerySeq', 'Iteration', function(x, max = FALSE) {
  getQuerySeq(x@hits, max = max)
})

## @return list<BStringSet>
setMethod('getHitSeq', 'Iteration', function(x, max = FALSE) {
  getHitSeq(x@hits, max = max)
})

## @return list<BStringSet>
setMethod('getMatch', 'Iteration', function(x, max = FALSE) {
  getMatch(x@hits, max = max)
})

## PercIdentity, PercPositive, PercGaps ####

setMethod('getPercIdentity', 'Iteration', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## getMaxPercIdentity() does NOT yield the PercIdentity of the hsp with
## the highest PercIdenty, but the PercIdentity of the hsp with the highest
## bitscore !!!
setMethod('getMaxPercIdentity', 'Iteration', function(x) {
  vapply(x@hits, getMaxPercIdentity, 0L)
})

setMethod('getMaxPercIdentity', 'IterationList', function(x) {
  lapply(x, getMaxPercIdentity)
})

setMethod('getPercPositive', 'Iteration', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

setMethod('getPercGaps', 'Iteration', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryCoverage, HitCoverage ####

setMethod('getQueryCoverage', 'Iteration', function(x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["query_len"]]
})

setMethod('getQueryCoverage', 'IterationList', function(x) {
  lapply(x, getQueryCoverage)
})

setMethod('getHitCoverage', 'Iteration', function(x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["hit_len"]]
})

setMethod('getHitCoverage', 'IterationList', function(x) {
  lapply(x, getHitCoverage)
})

# subsetting, Iteration, IterationList -----------------------------------

## Subset to HitList
setMethod("[", "Iteration", function(x, i, j, ..., drop) {
  if (missing(i)) x@hits[] else x@hits[i]
})

## Subset to Hit
setMethod("[[", "Iteration", function(x, i, j, ...) {
  x@hits[[i]]
})


setMethod("is.na", "Iteration", function(x) {
  is(x@hits, "character") && x@hits == "No hits found"
})

setMethod("[", "IterationList", function(x, i, j, ..., drop) {
  IterationList(compact(callNextMethod()))
})

setMethod("[[", "IterationList", function(x, i, j, ...) {
  callNextMethod()
})

setMethod("is.na", "IterationList", function(x) {
  vapply(x, is.na, FALSE, USE.NAMES = FALSE)
})

# show, Iteration, IterationList -----------------------------------------

.show_Iteration <- function(it) {
  offset <- nchar(getIterNum(it)) + 8
  indent <- blanks(3)
  desc_ <- linebreak(getQueryDef(it), indent = -offset, offset = offset)
  header <- sprintf("Query %s: %s", getIterNum(it), desc_)
  
  if (is(it@hits, "character")) {
    tail <- NULL
    showme <- sprintf("%s%s\n", indent, it@hits)
  } 
  else {
    width <- getOption("width") - 44 - nchar(indent)
    if (width < 16) {
      showme <- paste0(header, "\n\nNot enough space to display hits")
    }
    else {
      tail <- "\n"
      n <- getOption("showHits", default = 12)
      assert_that(is.numeric(n), length(n) == 1, n > 0)
      nhits <- length(it@hits)
      if (n >= nhits) {
        n <- nhits
      }
      else {
        tail <- paste0(" ... and ", nhits - n, " more hits.\n")
      }
      n <- seq_len(n)
      hits <- getHit(it, n = n)
      desc <- getPrimaryHitDef(hits)
      mBSc <- getMaxBitscore(hits)
      tBSc <- getTotalBitscore(hits)
      qCov <- getQueryCoverage(hits)*100
      eVal <- getEvalue(hits, max = TRUE)
      accn <- getAccession(hits)
      c1 <- c(format("Description", width = width), blanks(width), format(ellipsize(desc, width), width=width))
      c2 <- c(" Max ", "Score", formatC(mBSc, digits=3, width=5, format='fg'))
      c3 <- c(" Total", " Score", formatC(tBSc, digits=3, width=6, format='fg'))
      c4 <- c(" Query", "  Cov ", paste0(formatC(qCov, digits=0, width=5, format='d'), '%'))
      c5 <- c("   E  ", " Value", formatC(eVal, digits=0, width=6, format='e'))
      c6 <- c("Accession   ", blanks(12), format(accn, width = 12))
      showme <- sprintf("%s%s %s %s %s %s  %s", indent, c1, c2, c3, c4, c5, c6)
    }
  }
  
  cat(header, showme, tail, sep="\n")
}

setMethod('show', 'Iteration',
          function(object) {
            if (is(object@hits, "character")) {
              olen <- "no"
            } else {
              olen <- length(object@hits)
            }        
            cat(sprintf("An %s instance with %s hit%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            .show_Iteration(object)
          })

setMethod('show', 'IterationList',
          function(object) {
            olen <- length(object)
            cat(sprintf("An %s instance with %s iteration%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            op <- options("showHits" = getOption("showHits", default = 3L))
            x <- lapply(object, .show_Iteration)
            options(op)
          })
