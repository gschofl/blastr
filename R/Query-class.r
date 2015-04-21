#' @include all-generics.r defline.r Hsp-class.r Hit-class.r
NULL


# Query-class, QueryList-class -------------------------------------------


setClassUnion("HitListOrChar", members = c("HitList", "character"))

#' Class \code{"Query"}
#' 
#' An S4 class that holds data parsed from an NCBI BLAST XML \sQuote{iteration} 
#' element. Each \code{"Query"} instance contains data from one query.
#'
#' @slot iter_num The number of the queries (iterations); <\code{integer}>.
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
#' showClass("Query")
new_Query <- 
  setClass(Class = "Query",
           slots = c(iter_num  = "integer",
                     query_id  = "character",
                     query_def = "character",
                     query_len = "integer",
                     hits      = "HitListOrChar",
                     query_env = "environment"),
           prototype = prototype(query_env = new.env(parent = emptyenv()))
  )

#' Class \code{"QueryList"} of \code{"Query"} objects
#' 
#' An S4 class that holds a list of \code{\linkS4class{Query}} objects.
#'
#' @slot query_env Shared container for \code{query_id},
#'    \code{query_def}, \code{query_len}, and \code{hit_len};
#'    \code{"environment"}.
#' @slot .Data Inherited from the \code{\link{list}} class.
#' @seealso \code{\linkS4class{BlastReport}}
#' @keywords classes
#' @export
#' @examples 
#' showClass("QueryList")
setClass(Class    = "QueryList",
         contains = "list",
         validity = listclassValidator("QueryList", "Query")
)

## constructor
QueryList <- listclassConstructor("QueryList", "Query")

# Hsp, nhsps, HspNum ####

## @return list<Hsp|HspList>
setMethod("getHsp", "Query", function(x, i, drop = TRUE) {
  lapply(x@hits, getHsp, i = i, drop = drop)
})

#' @describeIn Query Returns the number of HSPs for each hit; <\code{integer}>. 
setMethod("nhsps", "Query", function(x) {
  vapply(x@hits, nhsps, 0L)
})

#' @describeIn QueryList Returns the number of HSPs for each hit for each query;
#' <\code{list}> of <\code{integer}s>.
setMethod('nhsps', 'QueryList', function(x) {
  lapply(x, nhsps)
})

## Hit, nhits, QueryNum, QueryID, QueryDef, QueryLen ####

#' @describeIn Query Return a \code{\linkS4class{Hit}} or \code{\linkS4class{HitList}}
setMethod("getHit", "Query", function(x, i, drop = TRUE, ...) {
  hit <- if (missing(i)) x@hits[] else x@hits[i]
  if (drop && length(hit) == 1) {
    hit[[1]]
  } else hit
})

#' @describeIn QueryList Returns a list of \code{\linkS4class{Hit}}s or \code{\linkS4class{HitList}}s
setMethod("getHit", "QueryList", function(x, i, drop = TRUE, ...) {
  f <- if (missing(i)) getHit else Partial(getHit, i = i)
  lapply(x, f, drop = drop)
})

#' @describeIn Query Returns the number of hits; <\code{numeric}>. 
setMethod("nhits", "Query", function(x) length(x@hits))

#' @describeIn QueryList Returns the numbers of hits; <\code{list}> of <\code{numeric}s>.
setMethod("nhits", "QueryList", function(x) {
  vapply(x, nhits, 0L)
})

#' @describeIn Query Returns the query (iteration) number; <\code{integer}>. 
setMethod("getQueryNum", "Query", function(x) x@iter_num)

#' @describeIn QueryList Returns the query (iteration) numbers; <\code{integer}>. 
setMethod("getQueryNum", "QueryList", function(x) {
  vapply(x, getQueryNum, 0L)
})

#' @describeIn Query Returns the query ID; <\code{character}>. 
setMethod("getQueryID", "Query", function(x) x@query_id)

#' @describeIn QueryList Returns the query IDs; <\code{character}>.
setMethod("getQueryID", "QueryList", function(x) {
  vapply(x, getQueryID, "")
})

#' @describeIn Query Returns the query definition; <\code{character}>. 
setMethod("getQueryDef", "Query", function(x) x@query_def)

#' @describeIn QueryList Returns the query definitions; <\code{character}>.
setMethod("getQueryDef", "QueryList", function(x) {
  vapply(x, getQueryDef, "")
})

#' @describeIn Query Returns the query length; <\code{integer}>. 
setMethod("getQueryLen", "Query", function(x) x@query_len)

#' @describeIn QueryList Returns the query lengths; <\code{integer}>.
setMethod("getQueryLen", "QueryList", function(x) {
  vapply(x, getQueryLen, 0L)
})

## HitNum, HitLen, Accession, GeneID ####

## @return vector<integer>
setMethod("getHitNum", "Query", function(x) {
  vapply(x@hits, getHitNum, 0L)
})

## @return list<vector<integer>>
setMethod("getHitNum", "QueryList", function(x) {
  lapply(x, getHitNum)
})

## @return vector<integer>
setMethod("getHitLen", "Query", function(x) {
  vapply(x@hits, getHitLen, 0L)
})

## @return list<vector<integer>> 
setMethod("getHitLen", "QueryList", function(x) {
  lapply(x, getHitLen)
})

## @return vector<character>
setMethod("getAccession", "Query", function(x) {
  vapply(x@hits, getAccession, "")
})

## @return list<vector<character>>
setMethod("getAccession", "QueryList", function(x) {
  lapply(x, getAccession)
})

## @return vector<character>
setMethod("getGeneID", "Query", function(x) {
  vapply(x@hits, getGeneID, "")
})

## @return list<vector<character>>
setMethod("getGeneID", "QueryList", function(x) {
  lapply(x, getGeneID)
})

## HitID, HitDef, Defline, PrimaryHitDef, PrimaryDefline ####

## @return list<matrix<character>> or list<vector<integer>>
setMethod("getHitID", "Query", function(x, db = 'any') {
  lapply(x@hits, getHitID, db = db)
})

## @return list<vector<character>>
setMethod("getHitDef", "Query", function(x) {
  lapply(x@hits, getHitDef)
})

## @return list<vector<character>>
setMethod("getDefline", "Query", function(x) {
  lapply(x@hits, getDefline)
})

## @return vector<character>
setMethod("getPrimaryHitDef", "Query", function(x) {
  vapply(x@hits, getPrimaryHitDef, "")
})
## @return list<vector<character>>
setMethod("getPrimaryHitDef", "QueryList", function(x) {
  lapply(x, getPrimaryHitDef)
})

## @return vector<character>
setMethod("getPrimaryDefline", "Query", function(x) {
  vapply(x@hits, getPrimaryDefline, "")
})

## @return list<vector<character>>
setMethod("getPrimaryDefline", "QueryList", function(x) {
  lapply(x, getPrimaryDefline)
})

## Bitscore, Score, Evalue ####

## list<vector<numeric>>
setMethod("getBitscore", "Query",
          function(x, max = FALSE, sum = FALSE) {
            getBitscore(x@hits, max = max, sum = sum)
          })

## @return vector<numeric>
setMethod('getMaxBitscore', 'Query', function(x) {
  getMaxBitscore(x@hits)
})

## @return list<vector<numeric>>
setMethod('getMaxBitscore', 'QueryList', function(x) {
  lapply(x, getMaxBitscore)
})

## @return vector<numeric>
setMethod('getTotalBitscore', 'Query', function(x) {
  getTotalBitscore(x@hits)
})

## @return list<vector<numeric>>
setMethod('getTotalBitscore', 'QueryList', function(x) {
  lapply(x, getTotalBitscore)
})

## @return list<vector<numeric>>
setMethod('getScore', 'Query', function(x, max = FALSE) {
  getScore(x@hits, max = max)
})

## @return list<vector<numeric>> 
setMethod('getEvalue', 'Query', function(x, max = FALSE) {
  getEvalue(x@hits, max = max)
})

## Identity, Positive, Gaps, AlignLen ####

## @return list<vector<integer>> 
setMethod('getIdentity', 'Query', function(x, max = FALSE) {
  getIdentity(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getPositive', 'Query', function(x, max = FALSE) {
  getPositive(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getGaps', 'Query', function(x, max = FALSE) {
  getGaps(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getAlignLen', 'Query', function(x, max = FALSE) {
  getAlignLen(x@hits, max = max)
})

## QueryFrom, HitFrom, QueryTo, HitTo ####

## @return list<vector<integer>> 
setMethod('getQueryFrom', 'Query', function(x, max = FALSE) {
  getQueryFrom(x@hits, max = max)
})

## @return list<vector<integer>> 
setMethod('getQueryTo', 'Query', function(x, max = FALSE) {
  getQueryTo(x@hits, max = max)
})

## @return list<vector<integer> 
setMethod('getHitFrom', 'Query', function(x, max = FALSE) {
  getHitFrom(x@hits, max = max)
})

## @return list<vector<integer> 
setMethod('getHitTo', 'Query', function(x, max = FALSE) {
  getHitTo(x@hits, max = max)
})

## QueryFrame, HitFrame ####

## @return list<vector<integer>>
setMethod('getQueryFrame', 'Query', function(x, max = FALSE) {
  getQueryFrame(x@hits, max = max)
})

## @return list<vector<integer>>
setMethod('getHitFrame', 'Query', function(x, max = FALSE) {
  getHitFrame(x@hits, max = max)
})

## QueryRange, HitRange ####

## @return list<IRanges> 
setMethod("getQueryRange", "Query", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getQueryRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
})

## @return list<IRanges> 
setMethod("getHitRange", "Query", function(x, max = FALSE) {
  ans <- IRangesList( lapply(x@hits, getHitRange, max = max) )
  if (max)
    IRanges::unlist(ans)
  else ans
}) 

## getQuerySeq, getHitSeq, getMatch ####

## @return list<BStringSet>
setMethod('getQuerySeq', 'Query', function(x, max = FALSE) {
  getQuerySeq(x@hits, max = max)
})

## @return list<BStringSet>
setMethod('getHitSeq', 'Query', function(x, max = FALSE) {
  getHitSeq(x@hits, max = max)
})

## @return list<BStringSet>
setMethod('getMatch', 'Query', function(x, max = FALSE) {
  getMatch(x@hits, max = max)
})

## PercIdentity, PercPositive, PercGaps ####

setMethod('getPercIdentity', 'Query', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercIdentity, max = max)
  if (max)
    unlist(ans)
  else ans
})

## getMaxPercIdentity() does NOT yield the PercIdentity of the hsp with
## the highest PercIdenty, but the PercIdentity of the hsp with the highest
## bitscore !!!
setMethod('getMaxPercIdentity', 'Query', function(x) {
  vapply(x@hits, getMaxPercIdentity, 0L)
})

setMethod('getMaxPercIdentity', 'QueryList', function(x) {
  lapply(x, getMaxPercIdentity)
})

setMethod('getPercPositive', 'Query', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercPositive, max = max)
  if (max)
    unlist(ans)
  else ans
})

setMethod('getPercGaps', 'Query', function(x, max = FALSE) {
  ans <- lapply(x@hits, getPercGaps, max = max)
  if (max)
    unlist(ans)
  else ans
})

## QueryCoverage, HitCoverage ####

setMethod('getQueryCoverage', 'Query', function(x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["query_len"]]
})

setMethod('getQueryCoverage', 'QueryList', function(x) {
  lapply(x, getQueryCoverage)
})

setMethod('getHitCoverage', 'Query', function(x) {
  x <- x@hits
  ans <- Map(function(frame, from, to) {
    sum(.range(frame, from, to, width=TRUE))
  }, frame=getQueryFrame(x), from=getQueryFrom(x), to=getQueryTo(x))
  unlist(ans)/x@query_env[["hit_len"]]
})

setMethod('getHitCoverage', 'QueryList', function(x) {
  lapply(x, getHitCoverage)
})

# subsetting, Query, QueryList -----------------------------------

## Subset to HitList
setMethod("[", "Query", function(x, i, j, ..., drop) {
  if (missing(i)) x@hits[] else x@hits[i]
})

## Subset to Hit
setMethod("[[", "Query", function(x, i, j, ...) {
  x@hits[[i]]
})


setMethod("is.na", "Query", function(x) {
  is(x@hits, "character") && x@hits == "No hits found"
})

setMethod("[", "QueryList", function(x, i, j, ..., drop) {
  QueryList(compact(callNextMethod()))
})

setMethod("[[", "QueryList", function(x, i, j, ...) {
  callNextMethod()
})

setMethod("is.na", "QueryList", function(x) {
  vapply(x, is.na, FALSE, USE.NAMES = FALSE)
})

# show, Query, QueryList -----------------------------------------

.show_Query <- function(it) {
  offset <- nchar(getQueryNum(it)) + 8
  indent <- blanks(3)
  desc_ <- linebreak(getQueryDef(it), indent = -offset, offset = offset)
  header <- sprintf("Query %s: %s", getQueryNum(it), desc_)
  
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

setMethod('show', 'Query',
          function(object) {
            if (is(object@hits, "character")) {
              olen <- "no"
            } else {
              olen <- length(object@hits)
            }        
            cat(sprintf("An %s instance with %s hit%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, '', 's')),
                sep="")
            .show_Query(object)
          })

setMethod('show', 'QueryList',
          function(object) {
            olen <- length(object)
            cat(sprintf("An %s instance with %s quer%s.\n",
                        sQuote(class(object)), olen, ifelse(olen == 1, 'y', 'ies')),
                sep="")
            op <- options("showHits" = getOption("showHits", default = 3L))
            x <- lapply(object, .show_Query)
            options(op)
          })
