#' @include blastReport-class.r
#' @importClassesFrom RSQLite SQLiteConnection
#' @importClassesFrom RSQLite SQLiteObject
#' @importClassesFrom RSQLite dbObjectId
#' @importClassesFrom DBI DBIConnection
#' @importClassesFrom DBI DBIObject
NULL

# blastReportDB-class ----------------------------------------------------


#' blastReportDB class
#' 
#' blastReportDB is an S4 class that represents a connection to an SQLite
#' database holding blast records organised in three tables:
#' query, hit, and hsp
#' 
#' @name blastReportDB-class
#' @rdname blastReportDB-class
#' @exportClass blastReportDB
setClass('blastReportDB', contains='SQLiteConnection')


setValidity('blastReportDB', function (object) {
  if (!all(c("hit","hsp","query") %in% dbListTables(object)))
    return("Table missing from 'blastReportDB'")
  if (!all(c("query_id","query_def","query_len") %in% dbListFields(object, "query")))
    return("Field missing from table 'query'")
  if (!all(c("query_id","hit_id","hit_num","gene_id","accession",
             "definition","length") %in% dbListFields(object, "hit")))
    return("Field missing from table 'query'")
  if (!all(c("query_id","hit_id","hsp_id","hsp_num","bit_score",
             "score","evalue","query_from","query_to","hit_from",
             "hit_to","query_frame","hit_frame","identity","positive",
             "gaps","align_len","qseq","hseq","midline") 
           %in% dbListFields(object, "hsp")))   
    return("Field missing from table 'query'")
  TRUE
})


# constructor, blastReportDB-class ####
blastReportDB <- function(dbPath = "~/local/workspace/sqlite3/sample64.db") {
  assert_that(is.readable(dbPath))
  con <- db_connect(dbPath)
  new("blastReportDB", con)
}

.getQueryDef <- getterConstructor('query_def', 'query', 'query_id')
#' @rdname QueryDef-methods
#' @aliases getQueryDef,blastReportFb-method
setMethod("getQueryDef", "blastReportDB", function (x, id) .getQueryDef(x, id))

.getQueryLen <- getterConstructor('query_len', 'query', 'query_id', as='integer')
#' @rdname QueryLen-methods
#' @aliases getQueryLen,blastReportDB-method
setMethod("getQueryLen", "blastReportDB", function (x, id) .getQueryLen(x, id))

.getHitID <- getterConstructor('hit_id', 'hit', 'query_id',as='integer')
#' @rdname HitID-methods
#' @aliases getHitID,blastReportDB-method
setMethod("getHitID", "blastReportDB", function (x, id)  .getHitID(x, id))

.getHitNum <- getterConstructor('hit_num', 'hit', 'query_id',as='integer')
#' @rdname HitNum-methods
#' @aliases getHitNum,blastReportDB-method
setMethod("getHitNum", "blastReportDB", function (x, id) .getHitNum(x, id))

.getHitLen <- getterConstructor('length', 'hit', 'query_id',as='integer')
#' @rdname HitLen-methods
#' @aliases getHitLen,blastReportDB-method
setMethod("getHitLen", "blastReportDB", function (x, id) .getHitLen(x, id))

.getAccession <- getterConstructor('accession', 'hit', 'query_id')
#' @rdname Accession-methods
#' @aliases getAccession,blastReportDB-method
setMethod("getAccession", "blastReportDB", function (x, id) .getAccession(x, id))

.getGeneID <- getterConstructor('gene_id', 'hit', 'query_id',as='integer')
#' @rdname GeneID-methods
#' @aliases getGeneID,blastReportDB-method
setMethod("getGeneID", "blastReportDB", function (x, id) .getGeneID(x, id))

.getHitDef <- getterConstructor('definition', 'hit', 'query_id')
#' @rdname HitDef-methods
#' @aliases getHitDef,blastReportDB-method
setMethod("getHitDef", "blastReportDB", function (x, id) .getHitDef(x, id))

.getHspHitID <- getterConstructor('hit_id', 'hsp', 'query_id', as='integer')
#' @rdname HspHitID-methods
#' @aliases getHspHitID,blastReportDB-method
setMethod("getHspHitID", "blastReportDB", function (x, id) .getHspHitID(x, id))

.getHspID <- getterConstructor('hsp_id', 'hsp', 'query_id', as='integer')
#' @rdname HspID-methods
#' @aliases getHspID,blastReportDB-method
setMethod("getHspID", "blastReportDB", function (x, id) .getHspID(x, id))

.getHspNum <- getterConstructor('hsp_num', 'hsp', 'query_id', as='integer')
#' @usage getHspNum(x, id)
#' @rdname HspNum-methods
#' @aliases getHspNum,blastReportDB-method
#' @docType methods
setMethod("getHspNum", "blastReportDB", function (x, id) .getHspNum(x, id))

.getBitscore <- getterConstructor('bit_score', 'hsp', 'query_id', as='numeric')
#' @rdname Bitscore-methods
#' @aliases getBitscore,blastReportDB-method
setMethod("getBitscore", "blastReportDB", function (x, id, max = FALSE, sum=FALSE) {
  if (max) {
    .getMaxBitscore(x,id)
  } else if (sum) {
    .getTotalBitscore(x,id)
  } else {
    .getBitscore(x,id)
  }
})
.getMaxBitscore <- getterConstructor('MAX(bit_score)','hsp','query_id', as='numeric')
#' @rdname Bitscore-methods
#' @aliases getMaxBitscore,blastReportDB-method
setMethod("getMaxBitscore", "blastReportDB", function (x, id){
  .getMaxBitscore(x, id)
})
.getTotalBitscore <- getterConstructor('SUM(bit_score)','hsp','query_id', as='numeric')
#' @rdname Bitscore-methods
#' @aliases getTotalBitscore,blastReportDB-method
setMethod("getTotalBitscore", "blastReportDB", function (x, id) {
  .getTotalBitscore(x, id)
})

.getScore <- getterConstructor('score', 'hsp', 'query_id', as='integer')
.getMaxScore <- getterConstructor('MAX(score)', 'hsp', 'query_id', as='integer')
#' @rdname Score-methods
#' @aliases getScore,blastReportDB-method
setMethod("getScore", "blastReportDB", function (x, id,max=FALSE) {
  if (max) {
    .getMaxScore(x,id)
  } else {
    .getScore(x, id)
  }
})

.getEvalue <- getterConstructor('evalue', 'hsp', 'query_id', as='numeric')
.getMinEvalue <- getterConstructor('MIN(evalue)', 'hsp', 'query_id', as='numeric')
#' @rdname Evalue-methods
#' @aliases getEvalue,blastReportDB-method
setMethod("getEvalue", "blastReportDB", function (x, id, max=FALSE) {
  if (max) {
    .getMinEvalue(x,id)
  } else {
    .getEvalue(x, id)
  }
})

.getQueryFrom <- getterConstructor('query_from', 'hsp', 'query_id', as='integer')
#' @rdname QueryFrom-methods
#' @aliases getQueryFrom,blastReportDB-method
setMethod("getQueryFrom", "blastReportDB", function (x, id) .getQueryFrom(x, id))

.getQueryTo <- getterConstructor('query_to', 'hsp', 'query_id', as='integer')
#' @rdname QueryTo-methods
#' @aliases getQueryTo,blastReportDB-method
setMethod("getQueryTo", "blastReportDB", function (x, id) .getQueryTo(x, id),)

#' @rdname QueryRange-methods
#' @aliases getQueryRange,blastReportDB-method
setMethod("getQueryRange", "blastReportDB", function (x, id) 'not implemented')

.getHitFrom <- getterConstructor('hit_from', 'hsp', 'query_id', as='integer')
#' @rdname HitFrom-methods
#' @aliases getHitFrom,blastReportDB-method
setMethod("getHitFrom", "blastReportDB", function (x, id) .getHitFrom(x, id))

.getHitTo <- getterConstructor('hit_to', 'hsp', 'query_id', as='integer')
#' @rdname HitTo-methods
#' @aliases getHitTo,blastReportDB-method
setMethod("getHitTo", "blastReportDB", function (x, id) .getHitTo(x, id))

#' @rdname HitRange-methods
#' @aliases getHitRange,blastReportDB-method
setMethod("getHitRange", "blastReportDB", function (x, id) 'not implemented')

.getQueryFrame <- getterConstructor('query_frame', 'hsp', 'query_id', as='integer')
#' @rdname QueryFrame-methods
#' @aliases getQueryFrame,blastReportDB-method
setMethod("getQueryFrame", "blastReportDB", function (x, id) .getQueryFrame(x, id))

.getHitFrame <- getterConstructor('hit_frame', 'hsp', 'query_id', as='integer')
#' @rdname HitFrame-methods
#' @aliases getHitFrame,blastReportDB-method
setMethod("getHitFrame", "blastReportDB", function (x, id) .getHitFrame(x, id))

.getIdentity <- getterConstructor('identity', 'hsp', 'query_id', as='integer')
.getMaxIdentity <- getterConstructor('identity','hsp','query_id',
                                             FUN='MAX',VAL='bit_score', as='integer')
#' @rdname Identity-methods
#' @aliases getIdentity,blastReportDB-method
setMethod("getIdentity", "blastReportDB", function (x, id, max= FALSE) {
  if (max) {
    .getMaxIdentity(x,id)
  } else {
    .getIdentity(x, id)
  }
})


.getPositive <- getterConstructor('positive', 'hsp', 'query_id', as='integer')
.getMaxPositive <- getterConstructor('positive', 'hsp', 'query_id',
                                     FUN='MAX', VAL='bit_score', as='integer')
#' @rdname Positive-methods
#' @aliases getPositive,blastReportDB-method
setMethod("getPositive", "blastReportDB", function (x, id, max=FALSE) {
  if (max) {
    unlist( .getMaxPositive(x, id) )
  } else {
    .getPositive(x, id)
  }
})

.getGaps <- getterConstructor('gaps', 'hsp', 'query_id', as='integer')
#' @rdname Gaps-methods
#' @aliases getGaps,blastReportDB-method
setMethod("getGaps", "blastReportDB", function (x, id) .getGaps(x, id))

.getAlignLen <- getterConstructor('align_len', 'hsp', 'query_id', as='integer')
#' @rdname AlignLen-methods
#' @aliases getAlignLen,blastReportDB-method
setMethod("getAlignLen", "blastReportDB", function (x, id) .getAlignLen(x, id))

.getQuerySeq <- getterConstructor('qseq', 'hsp', 'query_id')
#' @rdname QuerySeq-methods
#' @aliases getQuerySeq,blastReportDB-method
setMethod("getQuerySeq", "blastReportDB", function (x, id) .getQuerySeq(x, id))

.getHitSeq <- getterConstructor('hseq', 'hsp', 'query_id')
#' @rdname HitSeq-methods
#' @aliases getHitSeq,blastReportDB-method
setMethod("getHitSeq", "blastReportDB", function (x, id) .getHitSeq(x, id))

.getMatch <-getterConstructor('midline', 'hsp', 'query_id')
#' @rdname Match-methods
#' @aliases getMatch,blastReportDB-method
setMethod("getMatch", "blastReportDB", function (x, id) .getMatch(x, id))

# .getPercIdentity <- getterConstructor(SELECT='CAST(identity AS FLOAT)/CAST(align_len AS FLOAT)',
#                                       FROM=eval(substitute(paste('(SELECT identity,align_len FROM hsp WHERE query_id=', id, ')'))))
# 
# require(pryr)
# unenclose(.getPercIdentity)
# .getPercIdentity(con, 16, TRUE)
# 
# 
# #' @rdname PercIdentity-methods
# #' @aliases getPercIdentity,blastReportDB-method
# setMethod("getPercIdentity", "blastReportDB", function (x, id) {
#   .getPercIdentity(x,id)
#   #as.list( unlist(getIdentity(x, id))/unlist(getAlignLen(x, id)) )
# })


#' @rdname PercIdentity-methods
#' @aliases getMaxPercIdentity,blastReportDB-method
setMethod("getMaxPercIdentity", "blastReportDB", function (x, id) {
  
})
  
## Noch überarbeiten
.getFrameFromAndTo <- getterConstructor("frame, from, to", "hsp", "query_id")

#' @rdname QueryCoverage-methods
#' @aliases getQueryCoverage,blastReportDB-method
setMethod("getQueryCoverage", "blastReportDB", function (x, id) {
  "Not implemented"
})

#' @rdname HitCoverage-methods
#' @aliases getHitCoverage,blastReportDB-method
setMethod("getHitCoverage", "blastReportDB", function (x, id) {
  "Not implemented"
})

#' @aliases show,blastReportDB-method
#' @rdname show-methods
setMethod('show', 'blastReportDB',
          function (object) {
            n <- db_count(object, "query")
            showme <- sprintf('%s object with %s query rows',
                              sQuote(class(object)), n)
            cat(showme, sep="\n")
          })
