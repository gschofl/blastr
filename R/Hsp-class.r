#' @include defline.r
#' @importFrom IRanges IRanges
#' @importFrom IRanges IRangesList
#' @importFrom IRanges reduce
#' @importFrom IRanges width
#' @importFrom IRanges unlist
#' @importFrom Biostrings BStringSet
NULL

setOldClass("list")


# Hsp-class --------------------------------------------------------------


#' Contains information parsed from an XML NCBI BLAST hsp element.
#'
#' Information about each high-scoring pair (hsp) within a hit contained in a
#' single \sQuote{Hsp} element.
#'
#'  \describe{
#'    \item{\code{hsp_num}:}{The number of the hsp; \code{"integer"}.}
#'    \item{\code{score}:}{The BLAST score of the hsp; \code{"numeric"}.}
#'    \item{\code{bit_score}:}{The bit score of the hsp; \code{"numeric"}.}
#'    \item{\code{evalue}:}{The expect value; \code{"numeric"}.}
#'    \item{\code{identity}:}{Number of identities; \code{"integer"}.}
#'    \item{\code{positive}:}{Number of positives; \code{"integer"}.}
#'    \item{\code{gaps}:}{Number of gaps; \code{"integer"}.}
#'    \item{\code{query_from}:}{Start residue for query sequence; \code{"integer"}.}
#'    \item{\code{query_to}:}{End residue for query sequence; \code{"integer"}.}
#'    \item{\code{hit_from}:}{Start residue for hit sequence; \code{"integer"}.}
#'    \item{\code{hit_to}:}{End residue for hit sequence; \code{"integer"}.}
#'    \item{\code{query_frame}:}{\code{"integer"}.}
#'    \item{\code{hit_frame}:}{\code{"integer"}.}
#'    \item{\code{qseq}:}{Query sequence; \code{"XString"}.}
#'    \item{\code{hseq}:}{Hit sequence; \code{"XString"}.}
#'    \item{\code{match}:}{Match sequence/midline; \code{"XString"}.}
#'    \item{\code{query_env}:}{Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.}
#'  } 
#' 
#' @name Hsp-class
#' @rdname Hsp-class
#' @exportClass Hsp
new_Hsp <- 
  setClass("Hsp",
           slots = c(hsp_num = "integer", score = "numeric",        
                     bit_score = "numeric", evalue = "numeric",
                     identity = "integer", positive = "integer",
                     gaps = "integer", align_len = "integer",
                     query_from = "integer", query_to = "integer",
                     hit_from = "integer", hit_to = "integer",
                     query_frame = "integer", hit_frame = "integer",
                     qseq = "XString", hseq = "XString",
                     match = "XString", query_env = 'environment')
  )

#' @name HspList-class
#' @rdname Hsp-class
#' @exportClass HspList
setClass("HspList", representation(query_env = 'environment'), contains="list",
         validity=listclassValidator('HspList', 'Hsp'))

## constructor
HspList <- listclassConstructor('HspList', 'Hsp')


# getter, Hsp, HspList ---------------------------------------------------

# yields index of the hsp with the highest bitscore
#
# except for score, evalue, and (obviously) bitscore we select the best
# hsp based on bitscore.
bs.max <- function(x) {
  which.max(vapply(x, slot, 'bit_score', FUN.VALUE=numeric(1)))
}

## HspNum ####

## For Hsp methods 'max' is a dummy argument which we need because the
## call to getHsp() can yield either a 'Hsp' or a 'HspList' object, so
## methods for both classes need provide the same arguments
##
##  : Hsp -> integer
##  : HspList [max = FALSE] -> vector<integer>
##  : HspList [max = TRUE] -> integer

#' @rdname HspNum-methods
#' @aliases getHspNum,Hsp-method
setMethod('getHspNum', 'Hsp', function(x, max = FALSE) x@hsp_num)

#' @rdname HspNum-methods
#' @aliases getHspNum,HspList-method
setMethod('getHspNum', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hsp_num', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## BitScore, Score, Evalue ####

##  : Hsp -> numeric
##  : HspList  [max,sum = FALSE] -> vector<numeric>
##  : HspList  [max,sum = TRUE] -> numeric
#' @rdname Bitscore-methods
#' @aliases getBitscore,Hsp-method
setMethod('getBitscore', 'Hsp', function(x, max = FALSE, sum = FALSE) x@bit_score)

#' @rdname Bitscore-methods
#' @aliases getBitscore,HspList-method
setMethod('getBitscore', 'HspList', function(x, max = FALSE, sum = FALSE) {
  ans <- vapply(x, slot, 'bit_score', FUN.VALUE=numeric(1))
  if (length(ans) > 1) {
    if (max) {
      ans[which.max(ans)]
    } else if (sum) {
      sum(ans)
    } else {
      ans
    }
  } else ans
})

#' @rdname Bitscore-methods
#' @aliases getBitscore,Hsp-method
setMethod('getMaxBitscore', 'Hsp', function(x) x@bit_score)

#' @rdname Bitscore-methods
#' @aliases getBitscore,HspList-method
setMethod('getMaxBitscore', 'HspList', function(x) getBitscore(x, max = TRUE))

#' @rdname Bitscore-methods
#' @aliases getBitscore,Hsp-method
setMethod('getTotalBitscore', 'Hsp', function(x) x@bit_score)

#' @rdname Bitscore-methods
#' @aliases getBitscore,HspList-method
setMethod('getTotalBitscore', 'HspList', function(x) getBitscore(x, sum = TRUE))

##  : Hsp -> numeric
##  : HspList  [max = FALSE] -> vector<numeric>
##  : HspList [max = TRUE] -> numeric

#' @rdname Score-methods
#' @aliases getScore,Hsp-method
setMethod('getScore', 'Hsp', function(x, max = FALSE) x@score)

#' @rdname Score-methods
#' @aliases getScore,HspList-method
setMethod('getScore', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'score', FUN.VALUE=numeric(1))
  if (max && length(ans) > 1L) {
    ans[which.max(ans)]
  }
  else ans
})

#' @rdname Evalue-methods
#' @aliases getEvalue,Hsp-method
setMethod('getEvalue', 'Hsp', function(x, max = FALSE) x@evalue)

#' @rdname Evalue-methods
#' @aliases getEvalue,HspList-method
setMethod('getEvalue', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'evalue', FUN.VALUE=numeric(1))
  if (max && length(ans) > 1L) {
    ans[which.min(ans)]
  }
  else ans
})

## Identity, Positive, Gaps, AlignLen ####

##  : Hsp -> integer
##  : HspList  [max = FALSE] -> vector<integer>
##  : HspList [max = TRUE] -> integer

#' @rdname Identity-methods
#' @aliases getIdentity,Hsp-method
setMethod('getIdentity', 'Hsp', function(x, max = FALSE) x@identity)

#' @rdname Identity-methods
#' @aliases getIdentity,HspList-method
setMethod('getIdentity', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'identity', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname Positive-methods
#' @aliases getPositive,Hsp-method
setMethod('getPositive', 'Hsp', function(x, max = FALSE) x@positive)

#' @rdname Positive-methods
#' @aliases getPositive,HspList-method
setMethod('getPositive', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'positive', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname Gaps-methods
#' @aliases getGaps,Hsp-method
setMethod('getGaps', 'Hsp', function(x, max = FALSE) x@gaps)

#' @rdname Gaps-methods
#' @aliases getGaps,HspList-method
setMethod('getGaps', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'gaps', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname AlignLen-methods
#' @aliases getAlignLen,Hsp-method
setMethod('getAlignLen', 'Hsp', function(x, max = FALSE) x@align_len)

#' @rdname AlignLen-methods
#' @aliases getAlignLen,HspList-method
setMethod('getAlignLen', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'align_len', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

##  : Hsp -> integer
##  : HspList  [max = FALSE] -> vector<integer>
##  : HspList [max = TRUE] -> integer

#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,Hsp-method
setMethod('getQueryFrom', 'Hsp', function(x, max = FALSE) x@query_from)

#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,HspList-method
setMethod('getQueryFrom', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_from', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname QueryTo-methods
#' @aliases getQueryTo,Hsp-method
setMethod('getQueryTo', 'Hsp', function(x, max = FALSE) x@query_to)

#' @rdname QueryTo-methods
#' @aliases getQueryTo,HspList-method
setMethod('getQueryTo', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_to', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname HitFrom-methods
#' @aliases getHitFrom,Hsp-method
setMethod('getHitFrom', 'Hsp', function(x, max = FALSE) x@hit_from)

#' @rdname HitFrom-methods
#' @aliases getHitFrom,HspList-method
setMethod('getHitFrom', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hit_from', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname HitTo-methods
#' @aliases getHitTo,Hsp-method
setMethod('getHitTo', 'Hsp', function(x, max = FALSE) x@hit_to)

#' @rdname HitTo-methods
#' @aliases getHitTo,HspList-method
setMethod('getHitTo', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hit_to', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## QueryFrame, HitFrame ####

##  : Hsp -> integer
##  : HspList  [max = FALSE] -> vector<integer>
##  : HspList [max = TRUE] -> integer

#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,Hsp-method
setMethod('getQueryFrame', 'Hsp', function(x, max = FALSE) x@query_frame)

#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,HspList-method
setMethod('getQueryFrame', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_frame', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## @return vector<integer> 
#' @rdname HitFrame-methods
#' @aliases getHitFrame,Hsp-method
setMethod('getHitFrame', 'Hsp', function(x, max = FALSE) x@hit_frame)

#' @rdname HitFrame-methods
#' @aliases getHitFrame,HspList-method
setMethod('getHitFrame', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hit_frame', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## QueryRange, HitRange ####


.range <- function(frame, from, to, width = FALSE) {
  start <- ifelse(frame >= 0L, from, to)
  end <- ifelse(frame >= 0L, to, from)
  r <- IRanges(start, end)
  if (width) width(reduce(r)) else r
}

##  : Hsp -> IRanges (single range)
##  : HspList  [max = FALSE] -> IRanges (possibly multiple ranges)
##  : HspList [max = TRUE] -> IRanges (single range)

#' @rdname QueryRange-methods
#' @aliases getQueryRange,Hsp-method
setMethod("getQueryRange", "Hsp", function(x, max = FALSE) {
  .range(x@query_frame, x@query_from, x@query_to)
})

#' @rdname QueryRange-methods
#' @aliases getQueryRange,HspList-method
setMethod("getQueryRange", "HspList", function(x, max = FALSE) {
  if (max)
    x <- x[bs.max(x)]
  .range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x))
})

#' @rdname HitRange-methods
#' @aliases getHitRange,Hsp-method
setMethod("getHitRange", "Hsp", function(x, max = FALSE) {
  .range(x@hit_frame, x@hit_from, x@hit_to)
})

#' @rdname HitRange-methods
#' @aliases getHitRange,HspList-method
setMethod("getHitRange", "HspList", function(x, max = FALSE) {
  if (max)
    x <- x[bs.max(x)]
  .range(getHitFrame(x), getHitFrom(x), getHitTo(x))
}) 

## getQuerySeq, getHitSeq, getMatch ####

##  : Hsp -> BString
##  : HspList -> BStringSet

#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,Hsp-method
setMethod('getQuerySeq', 'Hsp', function(x, max = FALSE) x@qseq)

#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,HspList-method
setMethod('getQuerySeq', 'HspList', function(x, max = FALSE) {
  ans <- BStringSet(lapply(x, slot, 'qseq'))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname HitSeq-methods
#' @aliases getHitSeq,Hsp-method
setMethod('getHitSeq', 'Hsp', function(x, max = FALSE) x@hseq)

#' @rdname HitSeq-methods
#' @aliases getHitSeq,HspList-method
setMethod('getHitSeq', 'HspList', function(x, max = FALSE) {
  ans <- BStringSet(lapply(x, slot, 'hseq'))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname Match-methods
#' @aliases getMatch,Hsp-method
setMethod('getMatch', 'Hsp', function(x, max = FALSE) x@match)

#' @rdname Match-methods
#' @aliases getMatch,HspList-method
setMethod('getMatch', 'HspList', function(x, max = FALSE) {
  ans <- BStringSet(lapply(x, slot, 'match'))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)][[1L]]
  }
  else ans
})

## PercIdentity, PercPositive, PercGaps ####

##  : Hsp -> numeric
##  : HspList  [max = FALSE] -> vector<numeric>
##  : HspList  [max = TRUE] -> numeric

#' @rdname PercIdentity-methods
#' @aliases getPercIdentity,Hsp-method
setMethod('getPercIdentity', 'Hsp', function(x, max = FALSE) {
  x@identity/x@align_len
})

#' @rdname PercIdentity-methods
#' @aliases getPercIdentity,HspList-method
setMethod('getPercIdentity', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, function(x) x@identity/x@align_len, numeric(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,Hsp-method
setMethod('getMaxPercIdentity', 'Hsp', function(x) {
  x@identity/x@align_len
})

#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,HspList-method
setMethod('getMaxPercIdentity', 'HspList', function(x) {
  best <- bs.max(x)
  getMaxPercIdentity(x[[best]])
})

##  : Hsp -> numeric
##  : HspList [max = FALSE] -> vector<numeric>
##  : HspList [max = TRUE] -> numeric
#' @rdname PercPositive-methods
#' @aliases getPercPositive,Hsp-method
setMethod('getPercPositive', 'Hsp', function(x, max = FALSE) {
  x@positive/x@align_len
})

#' @rdname PercPositive-methods
#' @aliases getPercPositive,HspList-method
setMethod('getPercPositive', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, function(x) x@positive/x@align_len, numeric(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

#' @rdname PercGaps-methods
#' @aliases getPercGaps,Hsp-method
setMethod('getPercGaps', 'Hsp', function(x, max = FALSE) {
  x@gaps/x@align_len
})

#' @rdname PercGaps-methods
#' @aliases getPercGaps,HspList-method
setMethod('getPercGaps', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, function(x) x@gaps/x@align_len, numeric(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## QueryCoverage, HitCoverage ####

##
## (Alignment length - gaps)/query[hit] length
## in case of overlapping hsps their ranges are
## normalized and the total alignment length is
## calculated
##
##  : Hsp -> numeric
##  : HspList [max = FALSE] -> vector<numeric>
##  : HspList [max = TRUE] -> numeric

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,Hsp-method
setMethod('getQueryCoverage', 'Hsp', function(x) {
  ( x@align_len - x@gaps )/x@query_env[['query_len']]
})

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,HspList-method
setMethod('getQueryCoverage', 'HspList', function(x) {
  sum(.range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x), width=TRUE))/
    x@query_env[['query_len']]
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,Hsp-method
setMethod('getHitCoverage', 'Hsp', function(x) {
  ( x@align_len - x@gaps )/x@query_env[['hit_len']]
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,HspList-method
setMethod('getHitCoverage', 'HspList', function(x) {
  sum(.range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x), width=TRUE))/
    x@query_env[['hit_len']]
})


# subsetting, HspList ----------------------------------------------------

setMethod("[", "HspList",
          function(x, i, j, ..., drop) {
            query_env <- x@query_env
            HspList( callNextMethod(), query_env = query_env )
          })


setMethod("[[", "HspList",
          function(x, i, j, ...) {
            callNextMethod() 
          })

# show, Hsp, HspList -----------------------------------------------------


.show_hsp <- function(hsp, offset = 0, show_aln = TRUE) {
  o <- blanks(offset)
  fmt <- paste0("\n%sRange %s: %s to %s\n",
                "%sScore: %s bits(%s), Expect: %s,\n",
                "%sIdentities: %s/%s(%s%%), Positives: %s/%s(%s%%), Gaps: %s/%s(%s%%)\n")
  head <- sprintf(fmt,
                  o, getHspNum(hsp), getHitFrom(hsp), getHitTo(hsp),
                  o, round(getBitscore(hsp), 1L), getScore(hsp),
                  format(getEvalue(hsp), scientific=TRUE, digits=2),
                  o, getIdentity(hsp), getAlignLen(hsp), round(100*getPercIdentity(hsp), 0),
                  getPositive(hsp), getAlignLen(hsp), round(100*getPercPositive(hsp), 0),
                  getGaps(hsp), getAlignLen(hsp), round(100*getPercGaps(hsp), 0))
  aln <- ""
  
  if (show_aln) {
    qstart <- min(getQueryFrom(hsp), getQueryTo(hsp))
    qrev <- ifelse(getQueryFrom(hsp) > getQueryTo(hsp), TRUE,  FALSE)
    hstart <- min(getHitFrom(hsp), getHitTo(hsp))
    hrev <- ifelse(getHitFrom(hsp) > getHitTo(hsp), TRUE, FALSE)
    aln <- wrapAlignment(toString(getQuerySeq(hsp)), toString(getMatch(hsp)), toString(getHitSeq(hsp)),
                         prefix=c("Query", "", "Spjct"), start=c(qstart, NA, hstart),
                         reverse=c(qrev, FALSE, hrev))
    aln <- paste0(aln, '\n')
  }
  
  cat(head, aln, sep="")
}

#' @aliases show,Hsp-method
#' @rdname show-methods
setMethod("show", "Hsp",
          function(object) {
            cat(sprintf("A %s instance", sQuote(class(object))))
            .show_hsp(object, show_aln = getOption("showAlignment", default=TRUE))
          })

#' @aliases show,HspList-method
#' @rdname show-methods
setMethod("show", "HspList",
          function(object) {
            olen <- length(object)
            tail <- "\n"
            n <- getOption("showHsps", default = 8)
            assert_that(is.numeric(n), length(n) == 1, n > 0)
            nhsps <- length(object)
            if (n >= nhsps) {
              n <- nhsps
            } else {
              tail <- paste0(" ... and ", nhsps - n, " more hsps.\n")
            }
            object <- object[seq_len(n)] 
            cat(sprintf("A %s instance with %s hsp%s.",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            show_aln <- logical(length(object))
            show_aln[] <- getOption("showAlignment", default=FALSE)
            x <- Map(function(hsp, sa) .show_hsp(hsp, show_aln=sa),
                     hsp = object, sa = as.list(show_aln))
          })
