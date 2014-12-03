#' @include blastReport-class.r
NULL


# subsetting, BlastReportDB ------------------------------------------------


## Subset to IterationList
setMethod("[", "BlastReportDB", function(x, i, j, ..., drop) {
  assert_that(is.numeric(i) | is.logical(i))
  if (is.logical(i)) {
    i <- which(i)
  }
  db <- BlastReportDB(db_path=":memory:", verbose=FALSE)
  WHERE <- paste0("where query_id in (", paste(i, collapse=","), ")")
  db_bulk_insert(db, "query", db_query(x, paste("select query_id, query_def, query_len from query", WHERE)), ...)
  db_bulk_insert(db, "hit", db_query(x, paste("select * from hit", WHERE)), ...)
  db_bulk_insert(db, "hsp", db_query(x, paste("select * from hsp", WHERE)), ...)
  db
})

## Subset to Iteration
setMethod("[[", "BlastReportDB", function(x, i, j, ...) {
  assert_that(length(i) == 1)
  x[i]
})

setMethod("is.na", "BlastReportDB", function(x) {
  db_query(db, 'select query_id from query', 1) %ni% db_query(db, 'select distinct query_id from hit', 1)
})

# getters, BlastReportDB --------------------------------------------------

setMethod("getQuery", "BlastReportDB", function(x, id, ...) {
  where <- if (!missing(id)) {
    paste0(" where query_id in (", paste0(id, collapse=","), ")")
  } else '' 
  stmt <- paste0("select query_id, query_def, query_len from query", where)
  db_query(x, stmt)
})

.getQueryID <- simpleGetter("query_id", 'query', WHERE='query_id', as='integer')
.getMaxQueryID <- simpleGetter("max(query_id)", 'query', WHERE='query_id', as='integer')
.getMinQueryID <- simpleGetter("min(query_id)", 'query', WHERE='query_id', as='integer')
setMethod("getQueryID", "BlastReportDB", function(x, id, ...) {
  .getQueryID(x, id, ...)
})


.getQueryDef <- simpleGetter('query_def', 'query', WHERE='query_id', as='character')
setMethod("getQueryDef", "BlastReportDB", function(x, id, ...) {
  .getQueryDef(x, id, ...)
})

.getQueryLen <- simpleGetter('query_len', 'query', WHERE='query_id', as='integer')
setMethod("getQueryLen", "BlastReportDB", function(x, id, ...) {
  .getQueryLen(x, id, ...)
})

.getHitID <- getterConstructor('hit_id', 'hit', WHERE='query_id', as='integer')
setMethod("getHitID", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getHitID(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getHitNum <- getterConstructor('hit_num', 'hit', WHERE='query_id', as='integer')
setMethod("getHitNum", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getHitNum(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getHitLen <- getterConstructor('length', 'hit', WHERE='query_id', as='integer')
setMethod("getHitLen", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getHitLen(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getAccession <- getterConstructor('accession', 'hit', WHERE='query_id', as='character')
setMethod("getAccession", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getAccession(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getGeneID <- getterConstructor('gene_id', 'hit', WHERE='query_id', as='character')
setMethod("getGeneID", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getGeneID(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getHitDef <- getterConstructor('definition', 'hit', WHERE='query_id', as='character')
setMethod("getHitDef", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  res <- .getHitDef(x, id, ...)
  if (length(res)==1) res[[1]] else res
})

.getHspHitID <- getterConstructor('hit_id', 'hsp', WHERE='query_id', as='integer')
.getMaxHspHitID <- getterConstructor('hit_id', 'hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', as='integer')
setMethod("getHspHitID", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max)
    unlist(.getMaxHspHitID(x, id, ...))
  else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHspHitID(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHspID <- getterConstructor('hsp_id', 'hsp', WHERE='query_id', as='integer')
.getMaxHspID <- getterConstructor('hsp_id', 'hsp', WHERE='query_id', FUN='max',
                                  VAL='bit_score', as='integer')
setMethod("getHspID", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max)
    unlist(.getMaxHspID(x, id, ...))
  else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHspID(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHspNum <- getterConstructor('hsp_num', 'hsp', WHERE='query_id', as='integer')
.getMaxHspNum <- getterConstructor('hsp_num', 'hsp', WHERE='query_id', FUN='max',
                                   VAL='bit_score', as='integer')
setMethod("getHspNum", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max)
    unlist(.getMaxHspNum(x, id, ...))
  else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHspNum(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getBitscore <- getterConstructor('bit_score', 'hsp', WHERE='query_id', as='numeric')
setMethod("getBitscore", "BlastReportDB", function(x, id, max=FALSE, sum=FALSE, ...) {
  if (max) {
    getMaxBitscore(x, id, ...)
  } else if (sum) {
    getTotalBitscore(x, id, ...)
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getBitscore(x, id)
    if (length(res)==1) res[[1]] else res
  }
})

.getMaxBitscore <- getterConstructor('max(bit_score)', 'hsp', WHERE='query_id',
                                     GROUP_BY='query_id', as='numeric')
setMethod("getMaxBitscore", "BlastReportDB", function(x, id, ...){
  unlist(.getMaxBitscore(x, id, ...))
})

.getTotalBitscore <- getterConstructor('sum(bit_score)', 'hsp', WHERE='query_id',
                                       GROUP_BY='query_id', as='numeric')
setMethod("getTotalBitscore", "BlastReportDB", function(x, id, ...) {
  unlist(.getTotalBitscore(x, id, ...))
})

.getScore <- getterConstructor('score', 'hsp', WHERE='query_id', as='integer')
.getMaxScore <- getterConstructor('max(score)', 'hsp', WHERE='query_id',
                                  GROUP_BY='query_id', as='integer')
setMethod("getScore", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxScore(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getScore(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getEvalue <- getterConstructor('evalue', 'hsp', WHERE='query_id', as='numeric')
.getMinEvalue <- getterConstructor('min(evalue)', 'hsp', WHERE='query_id',
                                   GROUP_BY='query_id', as='numeric')
setMethod("getEvalue", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMinEvalue(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getEvalue(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getQueryFrom <- getterConstructor('query_from', 'hsp', WHERE='query_id', as='integer')
.getMaxQueryFrom <- getterConstructor('query_from', 'hsp', WHERE='query_id', FUN='max',
                                      VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getQueryFrom", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxQueryFrom(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getQueryFrom(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getQueryTo <- getterConstructor('query_to', 'hsp', WHERE='query_id', as='integer')
.getMaxQueryTo <- getterConstructor('query_to', 'hsp', WHERE='query_id', FUN='max',
                                    VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getQueryTo", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxQueryTo(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getQueryTo(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHitFrom <- getterConstructor('hit_from', 'hsp', WHERE='query_id', as='integer')
.getMaxHitFrom <- getterConstructor('hit_from', 'hsp', WHERE='query_id', FUN='max',
                                    VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getHitFrom", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxHitFrom(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHitFrom(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHitTo <- getterConstructor('hit_to', 'hsp', WHERE='query_id', as='integer')
.getMaxHitTo <- getterConstructor('hit_to', 'hsp', WHERE='query_id', FUN='MAX',
                                  VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getHitTo", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxHitTo(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHitTo(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

setMethod("getQueryRange", "BlastReportDB", function(x, id, max = FALSE, ...) {
  res <- lapply(id, Compose(IRangesList, .rangeDB), con = conn(x), type = 'query', max = max, ...)
  if (length(res)==1) res[[1]] else res
})

setMethod("getHitRange", "BlastReportDB", function(x, id, max = FALSE, ...) {
  res <- lapply(id, Compose(IRangesList, .rangeDB), con = conn(x), type = 'hit', max = max, ...)
  if (length(res)==1) res[[1]] else res
})

.getQueryFrame <- getterConstructor('query_frame', 'hsp', WHERE='query_id', as='integer')
.getMaxQueryFrame <- getterConstructor('query_frame', 'hsp', WHERE='query_id', FUN='max',
                                       VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getQueryFrame", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxQueryFrame(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getQueryFrame(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHitFrame <- getterConstructor('hit_frame', 'hsp', WHERE='query_id', as='integer')
.getMaxHitFrame <- getterConstructor('hit_frame', 'hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getHitFrame", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxHitFrame(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHitFrame(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getIdentity <- getterConstructor('identity', 'hsp', WHERE='query_id', as='integer')
.getMaxIdentity <- getterConstructor('identity', 'hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getIdentity", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxIdentity(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getIdentity(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getPositive <- getterConstructor('positive', 'hsp', WHERE='query_id', as='integer')
.getMaxPositive <- getterConstructor('positive', 'hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getPositive", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxPositive(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getPositive(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getGaps <- getterConstructor('gaps', 'hsp', WHERE='query_id', as='integer')
.getMaxGaps <- getterConstructor('gaps', 'hsp', WHERE='query_id', FUN='max',
                                 VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getGaps", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxGaps(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getGaps(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getAlignLen <- getterConstructor('align_len', 'hsp', WHERE='query_id', as='integer')
.getMaxAlignLen <- getterConstructor('align_len','hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', GROUP_BY='query_id', as='integer')
setMethod("getAlignLen", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxAlignLen(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getAlignLen(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getQuerySeq <- getterConstructor('qseq', 'hsp', WHERE='query_id')
.getMaxQuerySeq <- getterConstructor('qseq', 'hsp', WHERE='query_id', FUN='max',
                                     VAL='bit_score', GROUP_BY='query_id', as='character')
setMethod("getQuerySeq", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxQuerySeq(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getQuerySeq(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getHitSeq <- getterConstructor('hseq', 'hsp', WHERE='query_id')
.getMaxHitSeq <- getterConstructor('hseq', 'hsp', WHERE='query_id', FUN='max',
                                   VAL='bit_score', GROUP_BY='query_id', as='character')
setMethod("getHitSeq", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxHitSeq(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getHitSeq(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getMatch <-getterConstructor('midline', 'hsp', WHERE='query_id')
.getMaxMatch <- getterConstructor('midline', 'hsp', WHERE='query_id', FUN='max',
                                  VAL='bit_score', GROUP_BY='query_id', as='character')
setMethod("getMatch", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxMatch(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getMatch(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

.getPercIdentity <- getterConstructor(SELECT='CAST(identity AS FLOAT)/CAST(align_len AS FLOAT)',
                                      FROM='hsp', WHERE='query_id ', as="numeric")
.getMaxPercIdentity <- getterConstructor(SELECT='MAX(CAST(identity AS FLOAT)/CAST(align_len AS FLOAT))',
                                         FROM='hsp', WHERE='query_id', FUN='max',
                                         VAL='bit_score', GROUP_BY='query_id', as="numeric")
setMethod("getPercIdentity", "BlastReportDB", function(x, id, max=FALSE, ...) {
  if (max) {
    unlist(.getMaxPercIdentity(x, id, ...))
  } else {
    if (missing(id))
      id <- db_query(x, "select query_id from query", 1L, ...)
    res <- .getPercIdentity(x, id, ...)
    if (length(res)==1) res[[1]] else res
  }
})

setMethod("getMaxPercIdentity", "BlastReportDB", function(x, id, ...) {
  unlist(.getMaxPercIdentity(x, id, ...))
})

setMethod("getQueryCoverage", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  hitwidth = lapply(.rangeDB(conn(x), id, type='query', width=TRUE, ...), FUN=vapply, "sum", 0L)
  querylen = .getQueryLen(x, id)
  res <- .mapply(`/`, list(hitwidth, querylen), NULL)
  if (length(res)==1) res[[1]] else res
}) 

setMethod("getHitCoverage", "BlastReportDB", function(x, id, ...) {
  if (missing(id))
    id <- db_query(x, "select query_id from query", 1L, ...)
  hitwidth = lapply(.rangeDB(conn(x), id, type='hit', width=TRUE, ...), FUN=vapply, "sum", 0L)
  querylen = .getQueryLen(x, id, ...)
  res <- .mapply(`/`, list(hitwidth, querylen), NULL)
  if (length(res)==1) res[[1]] else res
})
