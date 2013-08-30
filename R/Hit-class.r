#' @include Hsp-class.r
NULL

# Hit-class --------------------------------------------------------------


#' Hit-class
#' 
#' A container for information parsed from an XML NCBI BLAST hit element.
#' Information about multiple Blast hits within a query is contained in a
#' \code{HitList}.
#'  
#'  \describe{
#'    \item{\code{hit_num}:}{The number of the hit;  \code{"integer"}.}
#'    \item{\code{hit_def}:}{Hit definition; \code{"DeflineSet"}.}
#'    \item{\code{hit_acc}:}{Accession number; \code{"character"}.}
#'    \item{\code{hit_len}:}{Length of hit; \code{"integer"}.}
#'    \item{\code{hsps}:}{List of HSPs; \code{"\linkS4class{HspList}"}.}
#'    \item{\code{query_env}:}{Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.}
#'  }
#' 
#' @name Hit-class
#' @rdname Hit-class
#' @exportClass Hit
setClass("Hit",
         slots = c(hit_num = "integer",
                   hit_def = "DeflineSet",
                   hit_acc = "character",  
                   hit_len = "integer",
                   hsps = "HspList",
                   query_env = 'environment'))

#' @name HitList-class
#' @rdname Hit-class
#' @exportClass HitList
setClass("HitList", representation(query_env = 'environment'), contains="list",
         validity=listclassValidator('HitList', 'Hit'))

## constructor
HitList <- listclassConstructor('HitList', 'Hit')

## Hsp, nhsps, HspNum ####

## @return Hsp|HspList
#' @rdname Hsp-methods
#' @aliases getHsp,Hit-method
setMethod("getHsp", "Hit",
          function (x, n = NULL, drop = TRUE) {
            hsp <- if (is.null(n)) x@hsps else x@hsps[n]
            if (drop && length(hsp) == 1)
              hsp[[1]]
            else
              hsp
          })

## @return list<Hsp|HspList>
#' @rdname Hsp-methods
#' @aliases getHsp,HitList-method
setMethod("getHsp", "HitList", function (x, n = NULL, drop = TRUE) {
  lapply(x, getHsp, n = n, drop = drop)
})

## @return numeric
#' @rdname nhsps-methods
#' @aliases nhsps,Hit-method
setMethod('nhsps', 'Hit', function (x) length(x@hsps))

## @return vector<numeric>
#' @rdname nhsps-methods
#' @aliases nhsps,HitList-method
setMethod('nhsps', 'HitList', function (x) {
  vapply(x, nhsps, numeric(1))
})

## @return integer
#' @rdname HspNum-methods
#' @aliases getHspNum,Hit-method
setMethod('getHspNum', 'Hit', function (x, max = FALSE) {
  getHspNum(x@hsps, max = max)
})

## @return vector<integer> 
#' @rdname HspNum-methods
#' @aliases getHspNum,HitList-method
setMethod('getHspNum', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getHspNum, max = max)
  if (max)
    unlist(ans)
  else ans
})

## HitNum, HitLen, Accession, GeneID ####

## @return integer
#' @rdname HitNum-methods
#' @aliases getHitNum,Hit-method
setMethod("getHitNum", "Hit", function (x) x@hit_num)

## @return vector<integer> 
#' @rdname HitNum-methods
#' @aliases getHitNum,HitList-method
setMethod("getHitNum", "HitList", function (x) {
  vapply(x, getHitNum, integer(1))
})

## @return integer
#' @rdname HitLen-methods
#' @aliases getHitLen,Hit-method
setMethod("getHitLen", "Hit", function (x) x@hit_len)

## @return vector<integer> 
#' @rdname HitLen-methods
#' @aliases getHitLen,HitList-method
setMethod("getHitLen", "HitList", function (x) {
  vapply(x, getHitLen, integer(1))
})

# @return character
#' @rdname Accession-methods
#' @aliases getAccession,Hit-method
setMethod("getAccession", "Hit", function (x) x@hit_acc)

## @return vector<character>
#' @rdname Accession-methods
#' @aliases getAccession,HitList-method
setMethod("getAccession", "HitList", function (x) {
  vapply(x, getAccession, character(1))
})

## @return character
#' @rdname GeneID-methods
#' @aliases getGeneID,Hit-method
setMethod("getGeneID", "Hit", function (x) {
  .getDeflineID(x@hit_def[[1L]], db = 'gi')
})

## @return vector<character> 
#' @rdname GeneID-methods
#' @aliases getGeneID,HitList-method
setMethod("getGeneID", "HitList", function (x) {
  vapply(x, getGeneID, character(1))
})

## HitID, HitDef, Defline, PrimaryHitDef, PrimaryDefline ####

## @return matrix<character> or vector<character>
#' @rdname HitID-methods
#' @aliases getHitID,Hit-method
setMethod("getHitID", "Hit", function (x, db = 'any') {
  .getDeflineID(x@hit_def, db = db)
})

## @return list<matrix<character>> or list<vector<integer>> 
#' @rdname HitID-methods
#' @aliases getHitID,HitList-method
setMethod("getHitID", "HitList", function (x, db = 'any') {
  lapply(x, getHitID, db = db)
})

## @return vector<character>
#' @rdname HitDef-methods
#' @aliases getHitDef,Hit-method
setMethod("getHitDef", "Hit", function (x) {
  .deflineDesc(x@hit_def)
})

## @return list<vector<intecharacterger>>
#' @rdname HitDef-methods
#' @aliases getHitDef,HitList-method
setMethod("getHitDef", "HitList", function (x) {
  lapply(x, getHitDef)
})


## @return vector<character>
#' @rdname Defline-methods
#' @aliases getDefline,Hit-method
setMethod("getDefline", "Hit", function (x) {
  paste0(.deflineID(x@hit_def),' ',.deflineDesc(x@hit_def))
})

## @return list<vector<intecharacterger>>
#' @rdname Defline-methods
#' @aliases getDefline,HitList-method
setMethod("getDefline", "HitList", function (x) {
  lapply(x, getDefline)
})

## @return character
#' @rdname HitDef-methods
#' @aliases getPrimaryHitDef,Hit-method
setMethod("getPrimaryHitDef", "Hit", function (x) {
  .deflineDesc(x@hit_def[[1L]])
})

## @return vector<intecharacterger>
#' @rdname HitDef-methods
#' @aliases getPrimaryHitDef,HitList-method
setMethod("getPrimaryHitDef", "HitList", function (x) {
  vapply(x, getPrimaryHitDef, character(1))
})


## @return character
#' @rdname Defline-methods
#' @aliases getPrimaryDefline,Hit-method
setMethod("getPrimaryDefline", "Hit", function (x) {
  paste0(.deflineID(x@hit_def[[1L]]),' ',.deflineDesc(x@hit_def[[1L]]))
})

## @return vector<character>
#' @rdname Defline-methods
#' @aliases getPrimaryDefline,HitList-method
setMethod("getPrimaryDefline", "HitList", function (x) {
  vapply(x, getPrimaryDefline, character(1))
})

## Bitscore, Score, Evalue ####

## @return vector<numeric>
#' @rdname Bitscore-methods
#' @aliases getBitscore,Hit-method
setMethod('getBitscore', 'Hit', function (x, max = FALSE, sum = FALSE) {
  getBitscore(x@hsps, max = max, sum = sum)
})

## @return list<vector<numeric>>
#' @rdname Bitscore-methods
#' @aliases getBitscore,HitList-method
setMethod('getBitscore', 'HitList', function (x, max = FALSE, sum = FALSE) {
  ans <- lapply(x, getBitscore, max = max, sum = sum)
  if (max || sum)
    unlist(ans)
  else ans
})

## @return numeric
#' @rdname Bitscore-methods
#' @aliases getMaxBitscore,Hit-method
setMethod('getMaxBitscore', 'Hit', function (x) {
  getMaxBitscore(x@hsps)
})

## @return vector<numeric>
#' @rdname Bitscore-methods
#' @aliases getMaxBitscore,HitList-method
setMethod('getMaxBitscore', 'HitList', function (x) {
  vapply(x, getMaxBitscore, numeric(1))
})

## @return numeric
#' @rdname Bitscore-methods
#' @aliases getTotalBitscore,Hit-method
setMethod('getTotalBitscore', 'Hit', function (x) {
  getTotalBitscore(x@hsps)
})

## @return vector<numeric>
#' @rdname Bitscore-methods
#' @aliases getTotalBitscore,HitList-method
setMethod('getTotalBitscore', 'HitList', function (x) {
  vapply(x, getTotalBitscore, numeric(1))
})

## @return vector<numeric>
#' @rdname Score-methods
#' @aliases getScore,Hit-method
setMethod('getScore', 'Hit', function (x, max = FALSE) getScore(x@hsps, max = max))

## @return list<vector<numeric>>
#' @rdname Score-methods
#' @aliases getScore,HitList-method
setMethod('getScore', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getScore, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<numeric> 
#' @rdname Evalue-methods
#' @aliases getEvalue,Hit-method
setMethod('getEvalue', 'Hit', function (x, max = FALSE) {
  getEvalue(x@hsps, max = max)
})

## @return list<vector<numeric>> 
#' @rdname Evalue-methods
#' @aliases getEvalue,HitList-method
setMethod('getEvalue', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getEvalue, max = max)
  if (max)
    unlist(ans)
  else ans
})

## Identity, Positive, Gaps, AlignLen ####

## @return <vector<integer> 
#' @rdname Identity-methods
#' @aliases getIdentity,Hit-method
setMethod('getIdentity', 'Hit', function (x, max = FALSE) getIdentity(x@hsps, max = max))

## @return list<vector<integer>> 
#' @rdname Identity-methods
#' @aliases getIdentity,HitList-method
setMethod('getIdentity', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
#' @rdname Positive-methods
#' @aliases getPositive,Hit-method
setMethod('getPositive', 'Hit', function (x, max = FALSE) getPositive(x@hsps, max = max))

## @return list<vector<integer>> 
#' @rdname Positive-methods
#' @aliases getPositive,HitList-method
setMethod('getPositive', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
#' @rdname Gaps-methods
#' @aliases getGaps,Hit-method
setMethod('getGaps', 'Hit', function (x, max = FALSE) getGaps(x@hsps, max = max))

## @return list<vector<integer>> 
#' @rdname Gaps-methods
#' @aliases getGaps,HitList-method
setMethod('getGaps', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
#' @rdname AlignLen-methods
#' @aliases getAlignLen,Hit-method
setMethod('getAlignLen', 'Hit', function (x, max = FALSE) getAlignLen(x@hsps, max = max))

## @return list<vector<integer>> 
#' @rdname AlignLen-methods
#' @aliases getAlignLen,HitList-method
setMethod('getAlignLen', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getAlignLen, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

## @return vector<integer> 
#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,Hit-method
setMethod('getQueryFrom', 'Hit', function (x, max = FALSE) {
  getQueryFrom(x@hsps, max = max)
})

## @return list<vector<integer>> 
#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,HitList-method
setMethod('getQueryFrom', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getQueryFrom, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
#' @rdname QueryTo-methods
#' @aliases getQueryTo,Hit-method
setMethod('getQueryTo', 'Hit', function (x, max = FALSE) {
  getQueryTo(x@hsps, max = max)
})

## @return list<vector<integer>> 
#' @rdname QueryTo-methods
#' @aliases getQueryTo,HitList-method
setMethod('getQueryTo', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getQueryTo, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
#' @rdname HitFrom-methods
#' @aliases getHitFrom,Hit-method
setMethod('getHitFrom', 'Hit', function (x, max = FALSE) {
  getHitFrom(x@hsps, max = max)
})

## @return list<vector<integer>> 
#' @rdname HitFrom-methods
#' @aliases getHitFrom,HitList-method
setMethod('getHitFrom', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getHitFrom, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
#' @rdname HitTo-methods
#' @aliases getHitTo,Hit-method
setMethod('getHitTo', 'Hit', function (x, max = FALSE) {
  getHitTo(x@hsps, max = max)
})

## @return list<vector<integer>> 
#' @rdname HitTo-methods
#' @aliases getHitTo,HitList-method
setMethod('getHitTo', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getHitTo, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryFrame, HitFrame ####

## @return vector<integer> 
#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,Hit-method
setMethod('getQueryFrame', 'Hit', function (x, max = FALSE) {
  getQueryFrame(x@hsps, max = max)
})

## @return list<vector<integer>>
#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,HitList-method
setMethod('getQueryFrame', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getQueryFrame, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
#' @rdname HitFrame-methods
#' @aliases getHitFrame,Hit-method
setMethod('getHitFrame', 'Hit', function (x, max = FALSE) {
  getHitFrame(x@hsps, max = max)
})

## @return list<vector<integer>>
#' @rdname HitFrame-methods
#' @aliases getHitFrame,HitList-method
setMethod('getHitFrame', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getHitFrame, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryRange, HitRange ####

## @return IRanges 
#' @rdname QueryRange-methods
#' @aliases getQueryRange,Hit-method
setMethod("getQueryRange", "Hit", function (x, max = FALSE) {
  if (max)
    x <- x[[bs.max(x@hsps)]]
  .range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x))
})

## @return IRangesList|IRanges 
#' @rdname QueryRange-methods
#' @aliases getQueryRange,Hit-method
setMethod("getQueryRange", "HitList", function (x, max = FALSE) {
  ans <- IRangesList( lapply(x, getQueryRange, max = max) )
  if (max)
    unlist(ans)
  else ans
}) 

## @return IRanges 
#' @rdname HitRange-methods
#' @aliases getHitRange,Hit-method
setMethod("getHitRange", "Hit", function (x, max = FALSE) {
  if (max)
    x <- x[[bs.max(x@hsps)]]
  .range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x))
}) 

## @return IRangesList|IRanges 
#' @rdname HitRange-methods
#' @aliases getHitRange,HitList-method
setMethod("getHitRange", "HitList", function (x, max = FALSE) {
  ans <- IRangesList( lapply(x, getHitRange, max = max) )
  if (max)
    unlist(ans)
  else ans
}) 

## getQuerySeq, getHitSeq, getMatch ####

## @return BStringSet
#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,Hit-method
setMethod('getQuerySeq', 'Hit', function (x, max = FALSE) {
  getQuerySeq(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,HitList-method
setMethod('getQuerySeq', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getQuerySeq, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## @return BStringSet
#' @rdname HitSeq-methods
#' @aliases getHitSeq,Hit-method
setMethod('getHitSeq', 'Hit', function (x, max = FALSE) {
  getHitSeq(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
#' @rdname HitSeq-methods
#' @aliases getHitSeq,HitList-method
setMethod('getHitSeq', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getHitSeq, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## @return BStringSet
#' @rdname Match-methods
#' @aliases getMatch,Hit-method
setMethod('getMatch', 'Hit', function (x, max = FALSE) {
  getMatch(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
#' @rdname Match-methods
#' @aliases getMatch,HitList-method
setMethod('getMatch', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getMatch, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## PercIdentity, PercPositive, PercGaps ####

#' @rdname PercIdentity-methods
#' @aliases getPercIdentity,Hit-method
setMethod('getPercIdentity', 'Hit', function (x, max = FALSE) {
  getIdentity(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

#' @rdname PercIdentity-methods
#' @aliases getPercIdentity,HitList-method
setMethod('getPercIdentity', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getPercIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## getMaxPercIdentity() does NOT yield the PercIdentity of the hsp with
## the highed PercIdenty, but the PercIdentity of the hsp with the highest
## bitscore !!!
#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,Hit-method
setMethod('getMaxPercIdentity', 'Hit', function (x) getMaxPercIdentity(x@hsps))

#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,HitList-method
setMethod('getMaxPercIdentity', 'HitList', function (x) {
  vapply(x, getMaxPercIdentity, numeric(1))
})

#' @rdname PercPositive-methods
#' @aliases getPercPositive,Hit-method
setMethod('getPercPositive', 'Hit', function (x, max = FALSE) {
  getPositive(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

#' @rdname PercPositive-methods
#' @aliases getPercPositive,HitList-method
setMethod('getPercPositive', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getPercPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

#' @rdname PercGaps-methods
#' @aliases getPercGaps,Hit-method
setMethod('getPercGaps', 'Hit', function (x, max = FALSE) {
  getGaps(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

#' @rdname PercGaps-methods
#' @aliases getPercGaps,HitList-method
setMethod('getPercGaps', 'HitList', function (x, max = FALSE) {
  ans <- lapply(x, getPercGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryCoverage, HitCoverage ####

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,Hit-method
setMethod('getQueryCoverage', 'Hit', function (x) {
  getQueryCoverage(x@hsps)
})

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,HitList-method
setMethod('getQueryCoverage', 'HitList', function (x) {
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["query_len"]]
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,Hit-method
setMethod('getHitCoverage', 'Hit', function (x) {
  getHitCoverage(x@hsps)
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,HitList-method
setMethod('getHitCoverage', 'HitList', function (x) {
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["hit_len"]]
})

## subsetting, Hit, HitList -----------------------------------------------


## Subset to HspList
setMethod("[", "Hit",
          function(x, i, j, ..., drop) {
            x@hsps[i]
          })

## Subset to Hsp
setMethod("[[", "Hit",
          function(x, i, j, ...) {
            x@hsps[[i]]
          })

setMethod("[", "HitList",
          function(x, i, j, ..., drop) {
            query_env = x@query_env
            HitList( callNextMethod(), query_env = query_env )
          })

setMethod("[[", "HitList",
          function(x, i, j, ...) {
            callNextMethod()
          })


# show, Hit, HitList -----------------------------------------------------


.show_Hit <- function (hit, show_hsps = TRUE) {
  offset <- nchar(getHitNum(hit)) + 6
  id_ <- .deflineID(hit@hit_def[[1L]])
  desc_ <- linebreak(.deflineDesc(hit@hit_def[[1L]]), indent=-offset, offset=offset)
  line1 <- sprintf("Hit %s: %s", getHitNum(hit), desc_)
  fmt2 <- "Id: %s, Length: %s, No. hsps: %s\n"
  s <- sprintf(fmt2, id_, getHitLen(hit), nhsps(hit))
  line2 <- linebreak(s, indent = offset, offset = offset)
  cat(line1, line2, sep="\n")
  if (show_hsps) {
    hsps <- hit@hsps
    show_aln <- logical(length(hsps))
    show_aln[] <- getOption("showAlignment", default=FALSE)
    x <- Map(function (hsp, sa) .show_hsp(hsp, show_aln=sa, offset=3),
             hsp = hsps, sa = as.list(show_aln))
  }
}

#' @aliases show,Hit-method
#' @rdname show-methods
setMethod("show", "Hit",
          function (object) {
            olen <- length(object@hsps)
            cat(sprintf("A %s instance with %s hsp%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            .show_Hit(object, show_hsps = TRUE)
          })

#' @aliases show,HitList-method
#' @rdname show-methods
setMethod("show", "HitList",
          function (object) {
            olen <- length(object)
            tail <- "\n"
            n <- getOption("showHits", default = 8)
            assert_that(is.numeric(n), length(n) == 1, n > 0)
            nhits <- length(object)
            if (n >= nhits) {
              n <- nhits
            } else {
              tail <- paste0(" ... and ", nhits - n, " more hits.\n")
            }
            object <- object[seq_len(n)] 
            cat(sprintf("A %s instance with %s hit%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            x <- lapply(object, .show_Hit, show_hsps = FALSE)
            cat(tail)
          })
