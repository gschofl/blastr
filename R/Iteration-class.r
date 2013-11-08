#' @include Hit-class.r
NULL

# Iteration-class, IterationList-class -----------------------------------


#' Iteration-class
#' 
#' An S4 class that holds data parsed from an NCBI BLAST XML \sQuote{Iteration} 
#' element. Each Iteration contains data about a query.
#'
#' @section Slots:
#' 
#'    \describe{
#'      \item{\code{iter_num}:}{The number of the iteration;  \code{"integer"}.}
#'      \item{\code{query_id}:}{Query identifier; \code{"character"}.}
#'      \item{\code{query_def}:}{Query definition; \code{"character"}.}
#'      \item{\code{query_len}:}{Query length; \code{"integer"}.}
#'      \item{\code{hits}:}{List of hits; \code{"\linkS4class{HitList}"}.}
#'      \item{\code{query_env}:}{Shared container for \code{query_id},
#'      \code{query_def}, \code{query_len}, and \code{hit_len};
#'      \code{"environment"}.}
#'    }
#' 
#' @name Iteration-class
#' @rdname Iteration-class
#' @exportClass Iteration
setClass("Iteration",
         slots = c(iter_num = "integer",
                   query_id = "character",
                   query_def = "character",
                   query_len = "integer",
                   hits = "HitList",
                   query_env = "environment"),
         prototype = prototype(query_env=new.env(parent=emptyenv())))


#' IterationList-class
#' 
#' A list of \linkS4class{Iteration} objects.
#'
#' @section Slots:
#'    \describe{
#'      \item{\code{query_env}:}{Shared container for \code{query_id},
#'          \code{query_def}, \code{query_len}, and \code{hit_len};
#'          \code{"environment"}.}
#'      \item{\code{.Data}:}{Inherited from the \code{\link{list}} class.}
#'    }
#'  
#' @name IterationList-class
#' @rdname IterationList-class
#' @exportClass IterationList
setClass("IterationList", contains="list",
         validity=listclassValidator('IterationList', 'Iteration'))

## constructor
IterationList <- listclassConstructor('IterationList', 'Iteration')

# Hsp, nhsps, HspNum ####

## @return list<Hsp|HspList>
#' @rdname Hsp-methods
#' @aliases getHsp,Iteration-method
setMethod("getHsp", "Iteration", function (x, n = NULL, drop = TRUE) {
  lapply(x@hits, getHsp, n = n, drop = drop)
})

## @return vector<numeric>
#' @rdname nhsps-methods
#' @aliases nhsps,Iteration-method
setMethod('nhsps', 'Iteration', function (x) {
  vapply(x@hits, nhsps, numeric(1))
})

## @return list<vector<numeric>>
#' @rdname nhsps-methods
#' @aliases nhsps,IterationList-method
setMethod('nhsps', 'IterationList', function (x) {
  lapply(x, nhsps)
})


## Hit, nhits, IterNum, QueryID, QueryDef, QueryLen ####

## @return Hit|HitList
#' @rdname Hit-methods
#' @aliases getHit,Iteration-method
setMethod("getHit", "Iteration",
          function (x, n = NULL, drop = TRUE, ...) {
            hit <- if (is.null(n)) x@hits else x@hits[n]
            if (drop && length(hit) == 1)
              hit[[1]]
            else
              hit
          })

## @return list<Hit|HitList>
#' @rdname Hit-methods
#' @aliases getHit,IterationList-method
setMethod("getHit", "IterationList",
          function (x, n = NULL, drop = TRUE, ...) {
            lapply(x, getHit, n = n, drop = drop)
          })

## @return numeric
#' @rdname nhits-methods
#' @aliases nhits,Iteration-method
setMethod('nhits', 'Iteration', function (x) length(x@hits))

## @return vector<numeric>
#' @rdname nhits-methods
#' @aliases nhits,IterationList-method
setMethod('nhits', 'IterationList', function (x) {
  vapply(x, nhits, numeric(1))
})

## @return integer
#' @rdname IterNum-methods
#' @aliases getIterNum,Iteration-method
setMethod("getIterNum", "Iteration", function (x) x@iter_num)

## @return vector<integer> 
#' @rdname IterNum-methods
#' @aliases getIterNum,IterationList-method
setMethod("getIterNum", "IterationList", function (x) {
  vapply(x, getIterNum, integer(1))
})

## @return character
#' @rdname QueryID-methods
#' @aliases getQueryID,Iteration-method
setMethod("getQueryID", "Iteration", function (x) x@query_id)

## @return vector<character> 
#' @rdname QueryID-methods
#' @aliases getQueryID,IterationList-method
setMethod("getQueryID", "IterationList", function (x) {
  vapply(x, getQueryID, character(1))
})

## @return character
#' @rdname QueryDef-methods
#' @aliases getQueryDef,Iteration-method
setMethod("getQueryDef", "Iteration", function (x) x@query_def)

## @return vector<character>
#' @rdname QueryDef-methods
#' @aliases getQueryDef,IterationList-method
setMethod("getQueryDef", "IterationList", function (x) {
  vapply(x, getQueryDef, character(1))
})

## @return integer
#' @rdname QueryLen-methods
#' @aliases getQueryLen,Iteration-method
setMethod("getQueryLen", "Iteration", function (x) x@query_len)

## @return vector<integer>
#' @rdname QueryLen-methods
#' @aliases getQueryLen,IterationList-method
setMethod("getQueryLen", "IterationList", function (x) {
  vapply(x, getQueryLen, integer(1))
})

## HitNum, HitLen, Accession, GeneID ####

## @return vector<integer>
#' @rdname HitNum-methods
#' @aliases HitNum,Iteration-method
setMethod("getHitNum", "Iteration", function (x) {
  vapply(x@hits, getHitNum, integer(1))
})

## @return list<vector<integer>>
#' @rdname HitNum-methods
#' @aliases getHitNum,IterationList-method
setMethod("getHitNum", "IterationList", function (x) {
  lapply(x, getHitNum)
})

## @return vector<integer>
#' @rdname HitLen-methods
#' @aliases getHitLen,Iteration-method
setMethod("getHitLen", "Iteration", function (x) {
  vapply(x@hits, getHitLen, integer(1))
})

## @return list<vector<integer>> 
#' @rdname HitLen-methods
#' @aliases getHitLen,IterationList-method
setMethod("getHitLen", "IterationList", function (x) {
  lapply(x, getHitLen)
})

## @return vector<character>
#' @rdname Accession-methods
#' @aliases getAccession,Iteration-method
setMethod("getAccession", "Iteration", function (x) {
  vapply(x@hits, getAccession, character(1))
})

## @return list<vector<character>>
#' @rdname Accession-methods
#' @aliases getAccession,IterationList-method
setMethod("getAccession", "IterationList", function (x) {
  lapply(x, getAccession)
})

## @return vector<character>
#' @rdname GeneID-methods
#' @aliases getGeneID,Iteration-method
setMethod("getGeneID", "Iteration", function (x) {
  vapply(x@hits, getGeneID, character(1))
})

## @return list<vector<character>>
#' @rdname GeneID-methods
#' @aliases getGeneID,IterationList-method
setMethod("getGeneID", "IterationList", function (x) {
  lapply(x, getGeneID)
})

## HitID, HitDef, Defline, PrimaryHitDef, PrimaryDefline ####

## @return list<matrix<character>> or list<vector<integer>>
#' @rdname HitID-methods
#' @aliases getHitID,Iteration-method
setMethod("getHitID", "Iteration", function (x, db = 'any') {
  lapply(x@hits, getHitID, db = db)
})

## @return list<vector<character>>
#' @rdname HitDef-methods
#' @aliases getHitDef,Iteration-method
setMethod("getHitDef", "Iteration", function (x) {
  lapply(x@hits, getHitDef)
})

## @return list<vector<character>>
#' @rdname Defline-methods
#' @aliases getDefline,Iteration-method
setMethod("getDefline", "Iteration", function (x) {
  lapply(x@hits, getDefline)
})

## @return vector<character>
#' @rdname HitDef-methods
#' @aliases getPrimaryHitDef,Iteration-method
setMethod("getPrimaryHitDef", "Iteration", function (x) {
  vapply(x@hits, getPrimaryHitDef, character(1))
})
## @return list<vector<character>>
#' @rdname HitDef-methods
#' @aliases getPrimaryHitDef,IterationList-method
setMethod("getPrimaryHitDef", "IterationList", function (x) {
  lapply(x, getPrimaryHitDef)
})

## @return vector<character>
#' @rdname Defline-methods
#' @aliases getPrimaryDefline,Iteration-method
setMethod("getPrimaryDefline", "Iteration", function (x) {
  vapply(x@hits, getPrimaryDefline, character(1))
})

## @return list<vector<character>>
#' @rdname Defline-methods
#' @aliases getPrimaryDefline,IterationList-method
setMethod("getPrimaryDefline", "IterationList", function (x) {
  lapply(x, getPrimaryDefline)
})

## Bitscore, Score, Evalue ####

## list<vector<numeric>>
#' @rdname Bitscore-methods
#' @aliases getBitscore,Iteration-method
setMethod("getBitscore", "Iteration",
          function (x, max = FALSE, sum = FALSE) {
            getBitscore(x@hits, max = max, sum = sum)
          })

## @return vector<numeric>
#' @rdname Bitscore-methods
#' @aliases getMaxBitscore,Iteration-method
setMethod('getMaxBitscore', 'Iteration', function (x) {
  getMaxBitscore(x@hits)
})

## @return list<vector<numeric>>
#' @rdname Bitscore-methods
#' @aliases getMaxBitscore,IterationList-method
setMethod('getMaxBitscore', 'IterationList', function (x) {
  lapply(x, getMaxBitscore)
})

## @return vector<numeric>
#' @rdname Bitscore-methods
#' @aliases getTotalBitscore,Iteration-method
setMethod('getTotalBitscore', 'Iteration', function (x) {
  getTotalBitscore(x@hits)
})
## @return list<vector<numeric>>
#' @rdname Bitscore-methods
#' @aliases getTotalBitscore,IterationList-method
setMethod('getTotalBitscore', 'IterationList', function (x) {
  lapply(x, getTotalBitscore)
})

## @return list<vector<numeric>>
#' @rdname Score-methods
#' @aliases getScore,Iteration-method
setMethod('getScore', 'Iteration', function (x, max = FALSE) {
  getScore(x@hits, max = max)
})

## @return list<vector<numeric>> 
#' @rdname Evalue-methods
#' @aliases getEvalue,Iteration-method
setMethod('getEvalue', 'Iteration', function (x, max = FALSE) {
  getEvalue(x@hits, max = max)
})

## Identity, Positive, Gaps, AlignLen ####

## @return list<vector<integer>> 
#' @rdname Identity-methods
#' @aliases getIdentity,Iteration-method
setMethod('getIdentity', 'Iteration', function (x, max = FALSE) {
  getIdentity(x@hits, max = max)
})

## @return list<vector<integer>> 
#' @rdname Identity-methods
#' @aliases getIdentity,Iteration-method
setMethod('getPositive', 'Iteration', function (x, max = FALSE) {
  getPositive(x@hits, max = max)
})

## @return list<vector<integer>> 
#' @rdname Gaps-methods
#' @aliases getGaps,Iteration-method
setMethod('getGaps', 'Iteration', function (x, max = FALSE) {
  getGaps(x@hits, max = max)
})

## @return list<vector<integer>> 
#' @rdname AlignLen-methods
#' @aliases getAlignLen,Hit-method
setMethod('getAlignLen', 'Iteration', function (x, max = FALSE) {
  getAlignLen(x@hits, max = max)
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

## @return list<vector<integer>> 
#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,Iteration-method
setMethod('getQueryFrom', 'Iteration', function (x, max = FALSE) {
  getQueryFrom(x@hits, max = max)
})

## @return list<vector<integer>> 
#' @rdname QueryTo-methods
#' @aliases getQueryTo,Iteration-method
setMethod('getQueryTo', 'Iteration', function (x, max = FALSE) {
  getQueryTo(x@hits, max = max)
})

## @return list<vector<integer> 
#' @rdname HitFrom-methods
#' @aliases getHitFrom,Iteration-method
setMethod('getHitFrom', 'Iteration', function (x, max = FALSE) {
  getHitFrom(x@hits, max = max)
})

## @return list<vector<integer> 
#' @rdname HitTo-methods
#' @aliases getHitTo,Iteration-method
setMethod('getHitTo', 'Iteration', function (x, max = FALSE) {
  getHitTo(x@hits, max = max)
})

## QueryFrame, HitFrame ####

## @return list<vector<integer>>
#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,Iteration-method
setMethod('getQueryFrame', 'Iteration', function (x, max = FALSE) {
  getQueryFrame(x@hits, max = max)
})

## @return list<vector<integer>>
#' @rdname HitFrame-methods
#' @aliases getHitFrame,Iteration-method
setMethod('getHitFrame', 'Iteration', function (x, max = FALSE) {
  getHitFrame(x@hits, max = max)
})

## QueryRange, HitRange ####

## @return list<IRanges> 
#' @rdname QueryRange-methods
#' @aliases getQueryRange,Iteration-method
setMethod("getQueryRange", "Iteration", function (x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getQueryRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
})

## @return list<IRanges> 
#' @rdname HitRange-methods
#' @aliases getHitRange,Iteration-method
setMethod("getHitRange", "Iteration", function (x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getHitRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
}) 

## getQuerySeq, getHitSeq, getMatch ####

## @return list<BStringSet>
#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,Iteration-method
setMethod('getQuerySeq', 'Iteration', function (x, max = FALSE) {
  getQuerySeq(x@hits, max = max)
})

## @return list<BStringSet>
#' @rdname HitSeq-methods
#' @aliases getHitSeq,Iteration-method
setMethod('getHitSeq', 'Iteration', function (x, max = FALSE) {
  getHitSeq(x@hits, max = max)
})

## @return list<BStringSet>
#' @rdname Match-methods
#' @aliases getMatch,Iteration-method
setMethod('getMatch', 'Iteration', function (x, max = FALSE) {
  getMatch(x@hits, max = max)
})

## PercIdentity, PercPositive, PercGaps ####

#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,Iteration-method
setMethod('getPercIdentity', 'Iteration', function (x, max = FALSE) {
  ans <- lapply(x@hits, getPercIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## getMaxPercIdentity() does NOT yield the PercIdentity of the hsp with
## the highest PercIdenty, but the PercIdentity of the hsp with the highest
## bitscore !!!
#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,Iteration-method
setMethod('getMaxPercIdentity', 'Iteration', function (x) {
  vapply(x@hits, getMaxPercIdentity, numeric(1))
})

#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,IterationList-method
setMethod('getMaxPercIdentity', 'IterationList', function (x) {
  lapply(x, getMaxPercIdentity)
})

#' @rdname PercPositive-methods
#' @aliases getPercPositive,Iteration-method
setMethod('getPercPositive', 'Iteration', function (x, max = FALSE) {
  ans <- lapply(x@hits, getPercPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

#' @rdname PercGaps-methods
#' @aliases getPercGaps,Iteration-method
setMethod('getPercGaps', 'Iteration', function (x, max = FALSE) {
  ans <- lapply(x@hits, getPercGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryCoverage, HitCoverage ####

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,Iteration-method
setMethod('getQueryCoverage', 'Iteration', function (x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["query_len"]]
})

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,IterationList-method
setMethod('getQueryCoverage', 'IterationList', function (x) {
  lapply(x, getQueryCoverage)
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,Iteration-method
setMethod('getHitCoverage', 'Iteration', function (x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["hit_len"]]
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,IterationList-method
setMethod('getHitCoverage', 'IterationList', function (x) {
  lapply(x, getHitCoverage)
})

# subsetting, Iteration, IterationList -----------------------------------

## Subset to HitList
setMethod("[", "Iteration",
          function(x, i, j, ..., drop) {
            x@hits[i]
          })

## Subset to Hit
setMethod("[[", "Iteration",
          function(x, i, j, ...) {
            x@hits[[i]]
          })

setMethod("[", "IterationList",
          function(x, i, j, ..., drop) {
            IterationList( callNextMethod() )
          })

setMethod("[[", "IterationList",
          function(x, i, j, ...) {
            callNextMethod()
          })

# show, Iteration, IterationList -----------------------------------------

.show_Iteration <- function (it) {
  if (is.string(it)) {
    heaf <- tail <- NULL
    showme <- it
  } else {
    offset <- nchar(getIterNum(it)) + 8
    indent <- blanks(3)
    desc_ <- linebreak(getQueryDef(it), indent = -offset, offset = offset)
    header <- sprintf("Query %s: %s", getIterNum(it), desc_)
    width <- getOption("width") - 44 - nchar(indent)
    if (width < 16) {
      showme <- paste0(header, "\n\nNot enough space to display hits")
    } else {
      tail <- "\n"
      n <- getOption("showHits", default = 12)
      assert_that(is.numeric(n), length(n) == 1, n > 0)
      nhits <- length(it@hits)
      if (n >= nhits) {
        n <- nhits
      } else {
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

#' @aliases show,Iteration-method
#' @rdname show-methods
setMethod('show', 'Iteration',
          function (object) {
            olen <- length(object@hits)
            cat(sprintf("An %s instance with %s hit%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            .show_Iteration(object)
          })

#' @aliases show,IterationList-method
#' @rdname show-methods
setMethod('show', 'IterationList',
          function (object) {
            olen <- length(object)
            cat(sprintf("An %s instance with %s iteration%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            op <- options("showHits" = getOption("showHits", default = 3L))
            x <- lapply(object, .show_Iteration)
            options(op)
          })
