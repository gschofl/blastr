#' @include all-generics.r defline.r
#' @importFrom IRanges IRanges IRangesList reduce width unlist
#' @importFrom Biostrings BStringSet
NULL

setOldClass("list")


# Hsp-class --------------------------------------------------------------


#' Contains information parsed from an XML NCBI BLAST hsp element.
#'
#' Information about each high-scoring pair (hsp) within a hit contained in a
#' single \sQuote{Hsp} element.
#'
#' @slot hsp_num The number of the hsp; \code{"integer"}.
#' @slot score The BLAST score of the hsp; \code{"numeric"}.
#' @slot bit_score The bit score of the hsp; \code{"numeric"}.
#' @slot evalue The expect value; \code{"numeric"}.
#' @slot identity Number of identities; \code{"integer"}.
#' @slot positive Number of positives; \code{"integer"}.
#' @slot gaps Number of gaps; \code{"integer"}.
#' @slot query_from Start residue for query sequence; \code{"integer"}.
#' @slot query_to End residue for query sequence; \code{"integer"}.
#' @slot hit_from Start residue for hit sequence; \code{"integer"}.
#' @slot hit_to End residue for hit sequence; \code{"integer"}.
#' @slot query_frame \code{"integer"}.
#' @slot hit_frame \code{"integer"}.
#' @slot qseq Query sequence; \code{"XString"}.
#' @slot hseq Hit sequence; \code{"XString"}.
#' @slot match Match sequence/midline; \code{"XString"}.
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @seealso \code{"\linkS4class{BlastReport}"} 
#' @export
new_Hsp <- 
  setClass(Class = "Hsp",
           slots = c(
             hsp_num = "integer", score = "numeric",        
             bit_score = "numeric", evalue = "numeric",
             identity = "integer", positive = "integer",
             gaps = "integer", align_len = "integer",
             query_from = "integer", query_to = "integer",
             hit_from = "integer", hit_to = "integer",
             query_frame = "integer", hit_frame = "integer",
             qseq = "XString", hseq = "XString",
             match = "XString", query_env = 'environment')
  )

#' HspList
#' 
#' A list of \code{"\linkS4class{Hsp}"} objects.
#'
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @slot .Data Inherited from the \code{\link{list}} class.
#' @seealso \code{"\linkS4class{BlastReport}"} 
#' @export
setClass("HspList", representation(query_env = 'environment'), contains = "list",
         validity = listclassValidator('HspList', 'Hsp'))

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
setMethod('getHspNum', 'Hsp', function(x, max = FALSE) x@hsp_num)

setMethod('getHspNum', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hsp_num', FUN.VALUE = integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  } else ans
})

## BitScore, Score, Evalue ####

##  : Hsp -> numeric
##  : HspList  [max,sum = FALSE] -> vector<numeric>
##  : HspList  [max,sum = TRUE] -> numeric
setMethod('getBitscore', 'Hsp', function(x, max = FALSE, sum = FALSE) x@bit_score)

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

setMethod('getMaxBitscore', 'Hsp', function(x) x@bit_score)

setMethod('getMaxBitscore', 'HspList', function(x) getBitscore(x, max = TRUE))

setMethod('getTotalBitscore', 'Hsp', function(x) x@bit_score)

setMethod('getTotalBitscore', 'HspList', function(x) getBitscore(x, sum = TRUE))

##  : Hsp -> numeric
##  : HspList  [max = FALSE] -> vector<numeric>
##  : HspList [max = TRUE] -> numeric
setMethod('getScore', 'Hsp', function(x, max = FALSE) x@score)

setMethod('getScore', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'score', FUN.VALUE=numeric(1))
  if (max && length(ans) > 1L) {
    ans[which.max(ans)]
  }
  else ans
})

setMethod('getEvalue', 'Hsp', function(x, max = FALSE) x@evalue)

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
setMethod('getIdentity', 'Hsp', function(x, max = FALSE) x@identity)

setMethod('getIdentity', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'identity', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getPositive', 'Hsp', function(x, max = FALSE) x@positive)

setMethod('getPositive', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'positive', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getGaps', 'Hsp', function(x, max = FALSE) x@gaps)

setMethod('getGaps', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'gaps', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getAlignLen', 'Hsp', function(x, max = FALSE) x@align_len)

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
setMethod('getQueryFrom', 'Hsp', function(x, max = FALSE) x@query_from)

setMethod('getQueryFrom', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_from', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getQueryTo', 'Hsp', function(x, max = FALSE) x@query_to)

setMethod('getQueryTo', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_to', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getHitFrom', 'Hsp', function(x, max = FALSE) x@hit_from)

setMethod('getHitFrom', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'hit_from', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getHitTo', 'Hsp', function(x, max = FALSE) x@hit_to)

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
setMethod('getQueryFrame', 'Hsp', function(x, max = FALSE) x@query_frame)

setMethod('getQueryFrame', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, slot, 'query_frame', FUN.VALUE=integer(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

## @return vector<integer> 
setMethod('getHitFrame', 'Hsp', function(x, max = FALSE) x@hit_frame)

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
setMethod("getQueryRange", "Hsp", function(x, max = FALSE) {
  .range(x@query_frame, x@query_from, x@query_to)
})

setMethod("getQueryRange", "HspList", function(x, max = FALSE) {
  if (max)
    x <- x[bs.max(x)]
  .range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x))
})

setMethod("getHitRange", "Hsp", function(x, max = FALSE) {
  .range(x@hit_frame, x@hit_from, x@hit_to)
})

setMethod("getHitRange", "HspList", function(x, max = FALSE) {
  if (max)
    x <- x[bs.max(x)]
  .range(getHitFrame(x), getHitFrom(x), getHitTo(x))
}) 

## getQuerySeq, getHitSeq, getMatch ####

##  : Hsp -> BString
##  : HspList -> BStringSet
setMethod('getQuerySeq', 'Hsp', function(x, max = FALSE) x@qseq)

setMethod('getQuerySeq', 'HspList', function(x, max = FALSE) {
  ans <- BStringSet(lapply(x, slot, 'qseq'))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getHitSeq', 'Hsp', function(x, max = FALSE) x@hseq)

setMethod('getHitSeq', 'HspList', function(x, max = FALSE) {
  ans <- BStringSet(lapply(x, slot, 'hseq'))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getMatch', 'Hsp', function(x, max = FALSE) x@match)

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
setMethod('getPercIdentity', 'Hsp', function(x, max = FALSE) {
  x@identity/x@align_len
})

setMethod('getPercIdentity', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, function(x) x@identity/x@align_len, numeric(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getMaxPercIdentity', 'Hsp', function(x) {
  x@identity/x@align_len
})

setMethod('getMaxPercIdentity', 'HspList', function(x) {
  best <- bs.max(x)
  getMaxPercIdentity(x[[best]])
})

##  : Hsp -> numeric
##  : HspList [max = FALSE] -> vector<numeric>
##  : HspList [max = TRUE] -> numeric
setMethod('getPercPositive', 'Hsp', function(x, max = FALSE) {
  x@positive/x@align_len
})

setMethod('getPercPositive', 'HspList', function(x, max = FALSE) {
  ans <- vapply(x, function(x) x@positive/x@align_len, numeric(1))
  if (max && length(ans) > 1L) {
    ans[bs.max(x)]
  }
  else ans
})

setMethod('getPercGaps', 'Hsp', function(x, max = FALSE) {
  x@gaps/x@align_len
})

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
setMethod('getQueryCoverage', 'Hsp', function(x) {
  ( x@align_len - x@gaps )/x@query_env[['query_len']]
})

setMethod('getQueryCoverage', 'HspList', function(x) {
  sum(.range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x), width=TRUE))/
    x@query_env[['query_len']]
})

setMethod('getHitCoverage', 'Hsp', function(x) {
  ( x@align_len - x@gaps )/x@query_env[['hit_len']]
})

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
    aln <- wrap_alignment(toString(getQuerySeq(hsp)), toString(getMatch(hsp)), toString(getHitSeq(hsp)),
                          prefix = c("Query", "", "Spjct"), start = c(qstart, NA, hstart),
                          reverse = c(qrev, FALSE, hrev))
    aln <- paste0(aln, '\n')
  }
  
  cat(head, aln, sep="")
}

setMethod("show", "Hsp",
          function(object) {
            cat(sprintf("A %s instance", sQuote(class(object))))
            .show_hsp(object, show_aln = getOption("showAlignment", default=TRUE))
          })

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
