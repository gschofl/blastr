#' @include utils.r
NULL

# The blastr API -------------------------------------------------------

##    Getters in blastReport
##      getHeader, getParams, getIteration, getHit, getIterNum,
##      getQueryID, getQueryDef, getQueryLen
##
##    Getters in Itertation, IterationList
##      getHit, getIterNum, getQueryID, getQueryDef, getQueryLen
##     
##    Getters in Iteration
##      getHsp, getHitNum, getHitAccession, getHitLen, getHitID
##
##    Getters in Hit, HitList
##      getHsp, getHitNum, getHitAccession, getHitLen, getHitID
##
##    Getters in Hsp, HspList, Hit
##      getBitscore, getMaxBitscore, getTotalBitscore, getScore
##      getEvalue, 
##
##    Getters/setters in blastReportDB
##
##
##    Getters/setters in blastTable
##
##
##    Show methods all classes
##      show
##
##    Subsetting Methods for ***List, ***Set classes
##      [
##


# Hsp-accessor generics -------------------------------------------------


#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the Hsp indices.
#' 
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return index of the Hsp with the highest
#' bit score.
#' @return An numeric vector or a list of numeric vectors.
#' @export
setGeneric("getHspNum", function (x, ...) standardGeneric("getHspNum"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the bit score, raw score, or expect values.
#' 
#' @usage getBitscore(x, max = FALSE, sum = FALSE)
#' @param x  A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return only the highest score.
#' @param sum If \code{TRUE} sum scores from all hsps of a Blast hit.
#' @return An numeric vector or a list of numeric vectors.
#' @rdname scores
#' @export
setGeneric("getBitscore", function (x, ...) standardGeneric("getBitscore"))

#' @usage getMaxBitscore(x)
#' @rdname scores
#' @export
setGeneric("getMaxBitscore", function (x, ...) standardGeneric("getMaxBitscore"))

#' @usage getTotalBitscore(x)
#' @rdname scores
#' @export
setGeneric("getTotalBitscore", function (x, ...) standardGeneric("getTotalBitscore"))

#' @usage getScore(x, max = FALSE)
#' @rdname scores
#' @export
setGeneric("getScore", function (x, ...) standardGeneric("getScore"))

#' @usage getEvalue(x, max = FALSE)
#' @rdname scores
#' @export
setGeneric("getEvalue", function (x, ...) standardGeneric("getEvalue"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the numbers of identities.
#' 
#' @usage getIdentity(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getIdentity", function (x, ...) standardGeneric("getIdentity"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the numbers of positives.
#' 
#' @usage getPositive(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getPositive", function (x, ...) standardGeneric("getPositive"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the numbers of gaps.
#' 
#' @usage getGaps(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getGaps", function (x, ...) standardGeneric("getGaps"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the alignment lengths.
#' 
#' @usage getIdentity(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getAlignLen", function (x, ...) standardGeneric("getAlignLen"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the position of the start, end, frame, or range of the residues of
#' query and hit sequences.
#' 
#' @usage getQueryFrom(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors.
#' @rdname coordinates
#' @export
setGeneric("getQueryFrom", function (x, ...) standardGeneric("getQueryFrom"))

#' @usage getQueryTo(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getQueryTo", function (x, ...) standardGeneric("getQueryTo"))

#' @usage getQueryFrame(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getQueryFrame", function (x, ...) standardGeneric("getQueryFrame"))

#' @usage getHitFrom(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getHitFrom", function (x, ...) standardGeneric("getHitFrom"))

#' @usage getHitTo(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getHitTo", function (x, ...) standardGeneric("getHitTo"))

#' @usage getHitFrame(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getHitFrame", function (x, ...) standardGeneric("getHitFrame"))

#' @usage getQueryRange(x, max = FALSE)
#' @return An \linkS4class{IRanges} object.
#' @rdname coordinates
#' @export
setGeneric("getQueryRange", function (x, ...) standardGeneric("getQueryRange"))

#' @usage getHitRange(x, max = FALSE)
#' @rdname coordinates
#' @export
setGeneric("getHitRange", function (x, ...) standardGeneric("getHitRange"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract the query or hit sequence.
#' 
#' @usage getQuerySeq(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An \linkS4class{BString} or \linkS4class{BStringSet} object.
#' @rdname sequence
#' @export
setGeneric("getQuerySeq", function (x, ...) standardGeneric("getQuerySeq"))

#' @usage getHitSeq(x, max = FALSE)
#' @rdname sequence
#' @export
setGeneric("getHitSeq", function (x, ...) standardGeneric("getHitSeq"))

#' @usage getMatch(x, max = FALSE)
#' @rdname sequence
#' @export
setGeneric("getMatch", function (x, ...) standardGeneric("getMatch"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract Percent Identity, Percent Positives or Percent Gaps.
#' 
#' @usage getPercIdentity(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getPercIdentity", function (x, ...) standardGeneric("getPercIdentity"))

#' @usage getMaxPercIdentity(x)
#' @rdname getPercIdentity
#' @export
setGeneric("getMaxPercIdentity", function (x, ...) standardGeneric("getMaxPercIdentity"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract Percent Positives.
#' 
#' @usage getPercPositive(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getPercPositive", function (x, ...) standardGeneric("getPercPositive"))

#' Access components of High-Scoring Pairs (Hsps)
#' 
#' Extract Percent Gaps.
#' 
#' @usage getPercGaps(x, max = FALSE)
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @param max If \code{TRUE} return data only for the hsp with the
#' highest bit score.
#' @return An numeric vector or a list of numeric vectors
#' @export
setGeneric("getPercGaps", function (x, ...) standardGeneric("getPercGaps"))

#' Extract the total Query Coverage.
#' 
#' Caluculated as (alignment length - gaps) / query length.
#' If hsps overlap, they are merged before caluclating the coverage. 
#' 
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @usage getQueryCoverage(x)
#' @export
setGeneric("getQueryCoverage", function (x, ...) standardGeneric("getQueryCoverage"))

#' Extract the total Hit Coverage.
#' 
#' Caluculated as (alignment length - gaps) / hit length.
#' If hsps overlap, they are merged before caluclating the coverage. 
#' 
#' @param x A \linkS4class{Hsp} or \linkS4class{HspList}.
#' @usage getHitCoverage(x)
#' @export
setGeneric("getHitCoverage", function (x, ...) standardGeneric("getHitCoverage"))


# Hit-accessor generics --------------------------------------------------


#' Access components of Blast Hits
#' 
#' Extract Hsp IDs
#' 
#' @usage getHspID(x)
#' @param x A \linkS4class{blastReportDB} connection.
#' @return An numeric vector.
#' @export
setGeneric("getHspID", function (x, ...) standardGeneric("getHspID"))

#' Access components of Blast Hits
#' 
#' Extract Hit IDs for Hsps. 
#' 
#' @usage getHspHitID(x)
#' @param x A \linkS4class{blastReportDB} connection.
#' @return An numeric vector.
#' @export
setGeneric("getHspHitID", function (x, ...) standardGeneric("getHspHitID"))

#' Access components of Blast Hits
#' 
#' Extract the number of Hsps for each Hit.
#' 
#' @usage nhsps(x)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @return An numeric vector.
#' @export
setGeneric("nhsps", function (x, ...) standardGeneric("nhsps"))

#' Access components of Blast Hits
#' 
#' Extract Hsps from a Hit.
#' 
#' @usage getHsp(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return An \linkS4class{Hsp} or \linkS4class{HspList}.
#' @export
setGeneric("getHsp", function (x, ...) standardGeneric("getHsp"))

#' Access components of Blast Hits
#' 
#' Extract the Hit indices.
#' 
#' @usage getHitNum(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return An \linkS4class{Hsp} or \linkS4class{HspList}.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getHitNum", function (x, ...) standardGeneric("getHitNum"))

#' Access components of Blast Hits
#' 
#' Extract the hit length.
#' 
#' @usage getHitLen(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getHitLen", function (x, ...) standardGeneric("getHitLen"))

#' Access components of Blast Hits
#' 
#' Extract hit IDs.
#' 
#' @usage getHitID(x, db = 'any', ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param db Database tag (e.g. 'gi', 'ref', 'gb', 'dbj', 'any')
#' @param ... Further arguments passed to methods.
#' @return An character matrix or character vector.
#' @export
setGeneric("getHitID", function (x, ...) standardGeneric("getHitID"))

#' Access components of Blast Hits
#' 
#' Extract the accession numbers of hits.
#' 
#' @usage getAccession(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return A character vector or list of character vectors.
#' @export
setGeneric("getAccession", function (x, ...) standardGeneric("getAccession"))

#' Access components of Blast Hits
#' 
#' Extract the GIs (gene identifiers) of hits.
#' 
#' @usage getGeneID(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return A character vector or list of character vectors.
#' @export
setGeneric("getGeneID", function (x, ...) standardGeneric("getGeneID"))

#' Access components of Blast Hits
#' 
#' Extract the (primary) hit definition.
#' 
#' @usage getHitDef(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return A character vector or list of character vectors.
#' @export
setGeneric("getHitDef", function (x, ...) standardGeneric("getHitDef"))

#' @usage getPrimaryHitDef(x, ...)
#' @rdname getHitDef
#' @export
setGeneric("getPrimaryHitDef", function (x, ...) standardGeneric("getPrimaryHitDef"))

#' Access components of Blast Hits
#' 
#' Extract the (primary) definition lines.
#' 
#' @usage getDefline(x, ...)
#' @param x A \linkS4class{Hit} or \linkS4class{HitList} object.
#' @param ... Further arguments passed to methods.
#' @return A character vector or list of character vectors.
#' @export
setGeneric("getDefline", function (x, ...) standardGeneric("getDefline"))

#' @usage getPrimaryDefline(x, ...)
#' @rdname getDefline
#' @export
setGeneric("getPrimaryDefline", function (x, ...) standardGeneric("getPrimaryDefline"))


# query-accessor generics ------------------------------------------------


#' Access components of Blast Iterations (Queries)
#' 
#' Extract a data frame with the fields:
#' query_id, query_def, query_len, n_hits
#' 
#' @usage getQuery(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{blastReportDB} object.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getQuery", function (x, ...) standardGeneric("getQuery"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract the Iteration (Query) numbers.
#' 
#' @usage getIterNum(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getIterNum", function (x, ...) standardGeneric("getIterNum"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract the query ID.
#' 
#' @usage getQueryID(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getQueryID", function (x, ...) standardGeneric("getQueryID"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract the query definition.
#' 
#' @usage getQueryDef(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @param ... Further arguments passed to methods.
#' @return A character vector or list of character vectors.
#' @export
setGeneric("getQueryDef", function (x, ...) standardGeneric("getQueryDef"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract the query length.
#' 
#' @usage getQueryLen(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector or list of numeric vectors.
#' @export
setGeneric("getQueryLen", function (x, ...) standardGeneric("getQueryLen"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract Blast hits.
#' 
#' @usage getHit(x, ...)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @param ... Further arguments passed to methods.
#' @return A \linkS4class{hit} or \linkS4class{HitList} object.
#' @export
setGeneric("getHit", function (x, ...) standardGeneric("getHit"))

#' Access components of Blast Iterations (Queries)
#' 
#' Extract the number of Hits for each Iteration/Query.
#' 
#' @usage nhits(x)
#' @param x A \linkS4class{blastReport} or \linkS4class{IterationList} object.
#' @return An numeric vector.
#' @export
setGeneric("nhits", function (x, ...) standardGeneric("nhits"))


# blastReport-accessor generics ------------------------------------------


#' Extract the header from a BLAST report
#' 
#' @usage getHeader(x, ...)
#' @param x A \code{\linkS4class{blastReport}} object.
#' @param ... Further arguments passed to methods.
#' @return A \linkS4class{blastHeader} object.
#' @export
setGeneric("getHeader", function (x, ...) standardGeneric("getHeader"))

#' Extract Iterations/Queries from a BLAST report
#' 
#' @usage getIteration(x, ...)
#' @param x A \code{\linkS4class{blastReport}} object.
#' @param ... Further arguments passed to methods.
#' @return A \linkS4class{Iteration} or \linkS4class{IterationList} object.
#' @export
setGeneric("getIteration", function (x, ...) standardGeneric("getIteration"))

#' Extract BLAST parameters and statistics from a report
#' 
#' @usage getParams(x, ...)
#' @param x A \code{\linkS4class{blastReport}} object.
#' @param ... Further arguments passed to methods.
#' @return A \linkS4class{blastParameters} object.
#' @export
setGeneric("getParams", function (x, ...) standardGeneric("getParams"))


# Defline, DeflineSet ----------------------------------------------------


#' @keywords internal
setGeneric(".deflineID", function (x, ...) standardGeneric(".deflineID"))

#' @keywords internal
setGeneric(".deflineDesc", function (x, ...) standardGeneric(".deflineDesc"))

#' @keywords internal
setGeneric(".getDeflineID", function (x, ...) standardGeneric(".getDeflineID"))

