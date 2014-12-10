#' @include all-generics.r Hsp-class.r defline.r
NULL

# Hit-class --------------------------------------------------------------


#' Class \code{"Hit"}
#' 
#' An S4 class that serves as a container for information parsed from an
#' NCBI BLAST XML hit element. Information about multiple BLAST hits within a 
#' query is contained in a \code{\linkS4class{HitList}}.
#'  
#' @slot hit_num The number of the hit; <\code{integer}>.
#' @slot hit_def Hit definition; <\code{\linkS4class{DeflineSet}}>.
#' @slot hit_acc Accession number; <\code{character}>.
#' @slot hit_len Length of hit; <\code{integer}>.
#' @slot hsps List of HSPs; <\code{\linkS4class{HspList}}>.
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    >\code{environment}>.
#' @seealso \code{\linkS4class{BlastReport}}
#' @export
#' @examples 
#' showClass("Hit")
new_Hit <- 
  setClass(Class = "Hit",
           slots = c(
             hit_num   = "integer",
             hit_def   = "DeflineSet",
             hit_acc   = "character",  
             hit_len   = "integer",
             hsps      = "HspList",
             query_env = 'environment'
           ))

#' Class \code{"HitList"} of \code{"Hit"} objects
#' 
#' An S4 class that holds a list of \code{\linkS4class{Hit}} objects.
#'
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @slot .Data Inherited from the \code{\link{list}} class.
#' @seealso \code{"\linkS4class{BlastReport}"}
#' @export
#' @examples 
#' showClass("HitList")
setClass(Class    = "HitList",
         slots    = c(query_env = 'environment'),
         contains = "list",
         validity = listclassValidator('HitList', 'Hit'))

## constructor
HitList <- listclassConstructor('HitList', 'Hit')

## Hsp, nhsps, HspNum ####

## @return Hsp|HspList
setMethod("getHsp", "Hit",
          function(x, n = NULL, drop = TRUE) {
            hsp <- if (is.null(n)) x@hsps else x@hsps[n]
            if (drop && length(hsp) == 1) hsp[[1]] else hsp
          })

## @return list<Hsp|HspList>
setMethod("getHsp", "HitList", function(x, n = NULL, drop = TRUE) {
  lapply(x, getHsp, n = n, drop = drop)
})

## @return numeric
setMethod('nhsps', 'Hit', function(x) length(x@hsps))

## @return vector<numeric>
setMethod('nhsps', 'HitList', function(x) {
  vapply(x, nhsps, numeric(1))
})

## @return integer
setMethod('getHspNum', 'Hit', function(x, max = FALSE) {
  getHspNum(x@hsps, max = max)
})

## @return vector<integer> 
setMethod('getHspNum', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getHspNum, max = max)
  if (max)
    unlist(ans)
  else ans
})

## HitNum, HitLen, Accession, GeneID ####

## @return integer
setMethod("getHitNum", "Hit", function(x) x@hit_num)

## @return vector<integer> 
setMethod("getHitNum", "HitList", function(x) {
  vapply(x, getHitNum, integer(1))
})

## @return integer
setMethod("getHitLen", "Hit", function(x) x@hit_len)

## @return vector<integer> 
setMethod("getHitLen", "HitList", function(x) {
  vapply(x, getHitLen, integer(1))
})

# @return character
setMethod("getAccession", "Hit", function(x) x@hit_acc)

## @return vector<character>
setMethod("getAccession", "HitList", function(x) {
  vapply(x, getAccession, character(1))
})

## @return character
setMethod("getGeneID", "Hit", function(x) {
  .getDeflineID(x@hit_def[[1L]], db = 'gi')
})

## @return vector<character> 
setMethod("getGeneID", "HitList", function(x) {
  vapply(x, getGeneID, character(1))
})

## HitID, HitDef, Defline, PrimaryHitDef, PrimaryDefline ####

## @return matrix<character> or vector<character>
setMethod("getHitID", "Hit", function(x, db = 'any') {
  .getDeflineID(x@hit_def, db = db)
})

## @return list<matrix<character>> or list<vector<integer>> 
setMethod("getHitID", "HitList", function(x, db = 'any') {
  lapply(x, getHitID, db = db)
})

## @return vector<character>
setMethod("getHitDef", "Hit", function(x) {
  .deflineDesc(x@hit_def)
})

## @return list<vector<intecharacterger>>
setMethod("getHitDef", "HitList", function(x) {
  lapply(x, getHitDef)
})


## @return vector<character>
setMethod("getDefline", "Hit", function(x) {
  paste0(.deflineID(x@hit_def),' ',.deflineDesc(x@hit_def))
})

## @return list<vector<intecharacterger>>
setMethod("getDefline", "HitList", function(x) {
  lapply(x, getDefline)
})

## @return character
setMethod("getPrimaryHitDef", "Hit", function(x) {
  .deflineDesc(x@hit_def[[1L]])
})

## @return vector<intecharacterger>
setMethod("getPrimaryHitDef", "HitList", function(x) {
  vapply(x, getPrimaryHitDef, character(1))
})

## @return character
setMethod("getPrimaryDefline", "Hit", function(x) {
  paste0(.deflineID(x@hit_def[[1L]]),' ',.deflineDesc(x@hit_def[[1L]]))
})

## @return vector<character>
setMethod("getPrimaryDefline", "HitList", function(x) {
  vapply(x, getPrimaryDefline, character(1))
})

## Bitscore, Score, Evalue ####

## @return vector<numeric>
setMethod('getBitscore', 'Hit', function(x, max = FALSE, sum = FALSE) {
  getBitscore(x@hsps, max = max, sum = sum)
})

## @return list<vector<numeric>>
setMethod('getBitscore', 'HitList', function(x, max = FALSE, sum = FALSE) {
  ans <- lapply(x, getBitscore, max = max, sum = sum)
  if (max || sum)
    unlist(ans)
  else ans
})

## @return numeric
setMethod('getMaxBitscore', 'Hit', function(x) {
  getMaxBitscore(x@hsps)
})

## @return vector<numeric>
setMethod('getMaxBitscore', 'HitList', function(x) {
  vapply(x, getMaxBitscore, numeric(1))
})

## @return numeric
setMethod('getTotalBitscore', 'Hit', function(x) {
  getTotalBitscore(x@hsps)
})

## @return vector<numeric>
setMethod('getTotalBitscore', 'HitList', function(x) {
  vapply(x, getTotalBitscore, numeric(1))
})

## @return vector<numeric>
setMethod('getScore', 'Hit', function(x, max = FALSE) getScore(x@hsps, max = max))

## @return list<vector<numeric>>
setMethod('getScore', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getScore, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<numeric> 
setMethod('getEvalue', 'Hit', function(x, max = FALSE) {
  getEvalue(x@hsps, max = max)
})

## @return list<vector<numeric>> 
setMethod('getEvalue', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getEvalue, max = max)
  if (max)
    unlist(ans)
  else ans
})

## Identity, Positive, Gaps, AlignLen ####

## @return <vector<integer> 
setMethod('getIdentity', 'Hit', function(x, max = FALSE) getIdentity(x@hsps, max = max))

## @return list<vector<integer>> 
setMethod('getIdentity', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
setMethod('getPositive', 'Hit', function(x, max = FALSE) getPositive(x@hsps, max = max))

## @return list<vector<integer>> 
setMethod('getPositive', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
setMethod('getGaps', 'Hit', function(x, max = FALSE) getGaps(x@hsps, max = max))

## @return list<vector<integer>> 
setMethod('getGaps', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return <vector<integer> 
setMethod('getAlignLen', 'Hit', function(x, max = FALSE) getAlignLen(x@hsps, max = max))

## @return list<vector<integer>> 
setMethod('getAlignLen', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getAlignLen, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

## @return vector<integer> 
setMethod('getQueryFrom', 'Hit', function(x, max = FALSE) {
  getQueryFrom(x@hsps, max = max)
})

## @return list<vector<integer>> 
setMethod('getQueryFrom', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getQueryFrom, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
setMethod('getQueryTo', 'Hit', function(x, max = FALSE) {
  getQueryTo(x@hsps, max = max)
})

## @return list<vector<integer>> 
setMethod('getQueryTo', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getQueryTo, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
setMethod('getHitFrom', 'Hit', function(x, max = FALSE) {
  getHitFrom(x@hsps, max = max)
})

## @return list<vector<integer>> 
setMethod('getHitFrom', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getHitFrom, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
setMethod('getHitTo', 'Hit', function(x, max = FALSE) {
  getHitTo(x@hsps, max = max)
})

## @return list<vector<integer>> 
setMethod('getHitTo', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getHitTo, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryFrame, HitFrame ####

## @return vector<integer> 
setMethod('getQueryFrame', 'Hit', function(x, max = FALSE) {
  getQueryFrame(x@hsps, max = max)
})

## @return list<vector<integer>>
setMethod('getQueryFrame', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getQueryFrame, max = max)
  if (max)
    unlist(ans)
  else ans
})

## @return vector<integer> 
setMethod('getHitFrame', 'Hit', function(x, max = FALSE) {
  getHitFrame(x@hsps, max = max)
})

## @return list<vector<integer>>
setMethod('getHitFrame', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getHitFrame, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryRange, HitRange ####

## @return IRanges 
setMethod("getQueryRange", "Hit", function(x, max = FALSE) {
  if (max)
    x <- x[[bs.max(x@hsps)]]
  .range(getQueryFrame(x), getQueryFrom(x), getQueryTo(x))
})

## @return IRangesList|IRanges 
setMethod("getQueryRange", "HitList", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x, getQueryRange, max = max) )
  if (max)
    unlist(ans)
  else ans
}) 

## @return IRanges 
setMethod("getHitRange", "Hit", function(x, max = FALSE) {
  if (max)
    x <- x[[bs.max(x@hsps)]]
  .range(getHitFrame(x), getHitFrom(x), getHitTo(x))
}) 

## @return IRangesList|IRanges 
setMethod("getHitRange", "HitList", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x, getHitRange, max = max) )
  if (max)
    unlist(ans)
  else ans
}) 

## getQuerySeq, getHitSeq, getMatch ####

## @return BStringSet
setMethod('getQuerySeq', 'Hit', function(x, max = FALSE) {
  getQuerySeq(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
setMethod('getQuerySeq', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getQuerySeq, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## @return BStringSet
setMethod('getHitSeq', 'Hit', function(x, max = FALSE) {
  getHitSeq(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
setMethod('getHitSeq', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getHitSeq, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## @return BStringSet
setMethod('getMatch', 'Hit', function(x, max = FALSE) {
  getMatch(x@hsps, max = max)
})

## @return list<BStringSet>|BStringSet
setMethod('getMatch', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getMatch, max = max)
  if (max)
    do.call('c', ans)
  else ans
})

## PercIdentity, PercPositive, PercGaps ####

setMethod('getPercIdentity', 'Hit', function(x, max = FALSE) {
  getIdentity(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

setMethod('getPercIdentity', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getPercIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## getMaxPercIdentity() does NOT yield the PercIdentity of the hsp with
## the highed PercIdenty, but the PercIdentity of the hsp with the highest
## bitscore !!!
setMethod('getMaxPercIdentity', 'Hit', function(x) getMaxPercIdentity(x@hsps))

setMethod('getMaxPercIdentity', 'HitList', function(x) {
  vapply(x, getMaxPercIdentity, numeric(1))
})

setMethod('getPercPositive', 'Hit', function(x, max = FALSE) {
  getPositive(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

setMethod('getPercPositive', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getPercPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

setMethod('getPercGaps', 'Hit', function(x, max = FALSE) {
  getGaps(x@hsps, max = max)/getAlignLen(x@hsps, max = max)
})

setMethod('getPercGaps', 'HitList', function(x, max = FALSE) {
  ans <- lapply(x, getPercGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryCoverage, HitCoverage ####

setMethod('getQueryCoverage', 'Hit', function(x) {
  getQueryCoverage(x@hsps)
})

setMethod('getQueryCoverage', 'HitList', function(x) {
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["query_len"]]
})

setMethod('getHitCoverage', 'Hit', function(x) {
  getHitCoverage(x@hsps)
})

setMethod('getHitCoverage', 'HitList', function(x) {
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["hit_len"]]
})

## subsetting, Hit, HitList -----------------------------------------------


## Subset to HspList
setMethod("[", "Hit", function(x, i, j, ..., drop) {
  x@hsps[i]
})

## Subset to Hsp
setMethod("[[", "Hit", function(x, i, j, ...) {
  x@hsps[[i]]
})

setMethod("[", "HitList", function(x, i, j, ..., drop) {
  query_env = x@query_env
  HitList( callNextMethod(), query_env = query_env )
})

setMethod("[[", "HitList", function(x, i, j, ...) {
  callNextMethod()
})


# show, Hit, HitList -----------------------------------------------------


.show_Hit <- function(hit, show_hsps = TRUE) {
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
    x <- Map(function(hsp, sa) .show_hsp(hsp, show_aln=sa, offset=3),
             hsp = hsps, sa = as.list(show_aln))
  }
}

setMethod("show", "Hit",
          function(object) {
            olen <- length(object@hsps)
            cat(sprintf("A %s instance with %s hsp%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            .show_Hit(object, show_hsps = TRUE)
          })

setMethod("show", "HitList",
          function(object) {
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
