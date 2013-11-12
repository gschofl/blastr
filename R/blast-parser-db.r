#' @include blastReportDB-class.r
#' @importFrom XML xmlEventParse
#' @importFrom XML xmlValue
#' @importFrom XML xmlStopParser
#' @importFrom XML xpathApply
NULL

BlastOutput.Iterations <- function(dbPath = NULL,
                                   max_hit = NULL,
                                   max_hsp = NULL,
                                   reset_at = 1000)
{
  if (!is.null(dbPath)) {
    con <- db_connect(dbPath)
    assert_that(con %has_tables% c("query", "hit", "hsp"))
    query_order <- dbListFields(con, "query")
    hit_order <- dbListFields(con, "hit")
    hsp_order <- dbListFields(con, "hsp")
  } else {
    query_order <- hit_order <- hsp_order <- TRUE
  }
  
  ## Environments to accumulate parsed data
  query <- new.env()
  hit <- new.env()
  hsp <- new.env()
  
  ## To keep track of hit_id and hsp_id
  init_hit <- 0
  init_hsp <- 0
  
  if (is.null(max_hit)) {
    max_hit <- 'last()'
  }
  
  if (is.null(max_hsp)) {
    max_hsp <- 'last()'
  }
  
  resetEnv <- function(...) {
    env <- list(...)
    lapply(env, function(e) {
      remove(list=ls(envir=e), envir=e, inherits=FALSE)
    })
  }
  
  accumulate <- function(tag, value, env) {
    assign(tag, c(env[[tag]], value), env) 
  }
  
  geneId <- function(def) {
    strsplitN(def, "|", 2L, fixed=TRUE)  
  }
  
  getLast <- function(tag, env, init_val) {
    val <- tryCatch(get(tag, envir=env, inherits=TRUE),
                    error=function(e) NULL)
    val[length(val)] %||% init_val
  }
  
  getHspRuns <- function(hsp_num) {
    x <- rle(diff(hsp_num))
    rval <- x$values
    rlen <- x$lengths
    x <- rlen[rval >= 0]
    ifelse(rval[rval >= 0] == 0, x, x + 1)
  }
  
  minimum <- function(a, b) {
    if (b == 'last()')
      b <- a
    ifelse(a < b, a, b)
  }
  
  `%.%` <- function(lval, rval) paste0(lval, rval)
  
  getQuery <- function() {
    ans <- as.list(query)[query_order]
    as.data.frame(ans, stringsAsFactors=FALSE)  
  }
  
  getHit <- function() {
    ans <- as.list(hit)[hit_order]
    as.data.frame(ans, stringsAsFactors=FALSE) 
  }
  
  getHsp <- function() {
    ans <- as.list(hsp)[hsp_order]
    as.data.frame(ans, stringsAsFactors=FALSE) 
  }
  
  db_load <- function(con) {
    assert_that( db_bulk_insert(con, "query", getQuery()) )
    assert_that( db_bulk_insert(con, "hit", getHit()) )
    assert_that( db_bulk_insert(con, "hsp", getHsp()) )
  }
  
  'Iteration_iter-num' <- function(ctxt, node) {
    .id <- as.integer(xmlValue(node))
    if (.id%%reset_at == 0) {
      print('Reset at '%.%.id)
      if (!is.null(dbPath)) {
        init_hit <<- getLast('hit_id', hit)
        init_hsp <<- getLast('hsp_id', hsp)
        db_load(con)
        resetEnv(query, hit, hsp)
      } else {
        xmlStopParser(ctxt)
      }
    }
    accumulate('query_id', .id, query)
  }
  class(`Iteration_iter-num`) <- "XMLParserContextFunction"
  
  'Iteration_query-def' <- function(node) {
    accumulate('query_def', xmlValue(node), query)
  }
  
  'Iteration_query-len' <- function(node) {
    accumulate('query_len', as.integer( xmlValue(node) ), query)
  }
  
  'Iteration_hits' <- function(node) {
    ## Accumulate over Hits
    xp <- paste0('/Iteration_hits/Hit[position() <= ', max_hit, ']/')
    
    hit_num <- xvalue(node,  xp%.%'Hit_num', as="integer")
    if (all(is.na(hit_num))) return(NULL)
    accumulate("hit_num",    hit_num, hit)
    
    hit_query_id <- rep(getLast('query_id', query), length(hit_num))
    accumulate("query_id",   hit_query_id, hit)
    
    hit_id <- getLast('hit_id', hit, init_hit) + hit_num
    accumulate("hit_id",     hit_id, hit)
    
    accumulate("gene_id",    xvalue(node, xp%.%'Hit_id', fun=geneId), hit)
    accumulate("accession",  xvalue(node, xp%.%'Hit_accession'), hit)
    accumulate("definition", xvalue(node, xp%.%'Hit_def'), hit)
    accumulate("length",     xvalue(node, xp%.%'Hit_len', as='integer'), hit)
    
    ## Accumulate over Hsps
    xp1 <- paste0('/Iteration_hits/Hit[position() <= ', max_hit, ']/Hit_hsps')
    xp2 <- paste0('/Hsp[position() <= ', max_hsp, ']/')
    
    hsp_num <- xvalue(node, xp1 %.% xp2 %.% 'Hsp_num', as="integer")
    #print("Hsp_num: "%.%paste0(hsp_num, collapse=","))
    accumulate("hsp_num",     hsp_num, hsp)
    
    hsp_len <- minimum(
      vapply(xpathApply(node, xp1), xsize, 'Hsp', FUN.VALUE=numeric(1)),
      max_hsp
    )
    #print("Hsp_len: "%.%paste0(hsp_len, collapse=","))
    
    hsp_query_id <- rep(getLast('query_id', query), length(hsp_num))
    #print("Hsp_query_id: "%.%paste0(hsp_query_id, collapse=","))
    accumulate("query_id",    hsp_query_id, hsp)
    
    hsp_hit_id <- unlist(Map(rep, hit_id, hsp_len))
    #print("Hsp_hit_id: "%.%paste0(hsp_hit_id, collapse=","))
    accumulate("hit_id",      hsp_hit_id, hsp)
    
    hsp_id <-  getLast('hsp_id', hsp, init_hsp) + seq_along(hsp_num)
    #print("Hsp_id: "%.%paste0(hsp_id, collapse=","))
    accumulate("hsp_id",      hsp_id, hsp)
    
    xp <- xp1 %.% xp2
    accumulate('bit_score',   xvalue(node, xp%.%'Hsp_bit-score', as='numeric'), hsp)
    accumulate('score',       xvalue(node, xp%.%'Hsp_score', as='integer'), hsp)
    accumulate('evalue',      xvalue(node, xp%.%'Hsp_evalue', as='numeric'), hsp)
    
    accumulate('query_from',  xvalue(node, xp%.%'Hsp_query-from', as='integer'), hsp)
    accumulate('query_to',    xvalue(node, xp%.%'Hsp_query-to', as='integer'), hsp)
    accumulate('hit_from',    xvalue(node, xp%.%'Hsp_hit-from', as='integer'), hsp)
    accumulate('hit_to',      xvalue(node, xp%.%'Hsp_hit-to', as='integer'), hsp)
    
    accumulate('query_frame', xvalue(node, xp%.%'Hsp_query-frame', as='integer'), hsp)
    accumulate('hit_frame',   xvalue(node, xp%.%'Hsp_hit-frame', as='integer'), hsp)
    accumulate('identity',    xvalue(node, xp%.%'Hsp_identity', as='integer'), hsp)
    accumulate('positive',    xvalue(node, xp%.%'Hsp_positive', as='integer'), hsp)
    accumulate('gaps',        xvalue(node, xp%.%'Hsp_gaps', as='integer'), hsp)
    accumulate('align_len',   xvalue(node, xp%.%'Hsp_align-len', as='integer'), hsp)
    
    accumulate('qseq',        xvalue(node, xp%.%'Hsp_qseq'), hsp)
    accumulate('hseq',        xvalue(node, xp%.%'Hsp_hseq'), hsp)
    accumulate('midline',     xvalue(node, xp%.%'Hsp_midline'), hsp)
  }
  
  list('Iteration_iter-num'=`Iteration_iter-num`,
       'Iteration_query-def'=`Iteration_query-def`,
       'Iteration_query-len'=`Iteration_query-len`,
       'Iteration_hits'=Iteration_hits,
       getQuery=getQuery,
       getHit=getHit,
       getHsp=getHsp)
}


blast_db.sql <- '
CREATE TABLE query(
        query_id      INTEGER,
        query_def     TEXT,
        query_len     INTEGER,
        PRIMARY KEY (query_id)
);
CREATE TABLE hit(
        query_id      INTEGER,
        hit_id        INTEGER,
        hit_num       INTEGER,
        gene_id       TEXT,
        accession     TEXT,
        definition    TEXT,
        length        INTEGER,
        PRIMARY KEY (hit_id),
        FOREIGN KEY (query_id) REFERENCES query (query_id)
);
CREATE TABLE hsp(
        query_id      INTEGER,
        hit_id        INTEGER,
        hsp_id        INTEGER,
        hsp_num       INTEGER,
        bit_score     FLOAT,
        score         INTEGER,
        evalue        FLOAT,
        query_from    INTEGER,
        query_to      INTEGER,
        hit_from      INTEGER,
        hit_to        INTEGER,
        query_frame   INTEGER,
        hit_frame     INTEGER,
        identity      INTEGER,
        positive      INTEGER,
        gaps          INTEGER,
        align_len     INTEGER,
        qseq          TEXT,
        hseq          TEXT,
        midline       TEXT,
        PRIMARY KEY (hsp_id),
        FOREIGN KEY (hit_id) REFERENCES hit (hit_id)
        FOREIGN KEY (query_id) REFERENCES query (query_id)
);
CREATE INDEX Fquery ON query (query_id);
CREATE INDEX Fhit ON hit (hit_id);
CREATE INDEX Fhit_query ON hit (query_id);
CREATE INDEX Fhit_hit_query ON hit (query_id, hit_id);
CREATE INDEX Fhsp ON hsp (hsp_id);
CREATE INDEX Fhsp_query ON hsp (query_id);
CREATE INDEX Fhsp_hit ON hsp (hit_id);
CREATE INDEX Fhsp_hit_query ON hsp (query_id, hit_id, hsp_id);
'

#' Load SQL schema for blast_db
#' 
#' @keywords internal
#' @export
.blastr_blastdb_sql <- function() {
  blast_db.sql
}


#' Parse NCBI BLAST XML files into \linkS4class{blastReportDB} objects.
#' 
#' Create (or connect to) a  blastReport SQLite database.
#' 
#' @param blastfile an XML BLAST Report.
#' @param db_path Path to an blastReport SQLite database.
#' @param max_hit How many hits should be parsed (default: all available)
#' @param max_hsp How many hsps should be parsed (default: all available)
#' @param reset_at After how many iterations should the parser dump the
#' data into the db before continuing.
#'
#' @return A \code{\linkS4class{blastReportDB}} object.
#' @rdname blastReportDB
#' @export
blastReportDB <- function(blastfile, db_path = "blast.db", max_hit = NULL,
                          max_hsp = NULL, reset_at = 1000)
{
  assert_that(is.readable(blastfile), has_extension(blastfile, 'xml'))
  con <- db_create(db_path, blast_db.sql, overwrite=TRUE)
  handler <- BlastOutput.Iterations(db_path, max_hit, max_hsp, reset_at)
  out <- xmlEventParse(blastfile, list(), branches=handler)
  ## load final part into db
  assert_that(db_bulk_insert(con, "query", handler$getQuery()))
  assert_that(db_bulk_insert(con, "hit", handler$getHit()))
  assert_that(db_bulk_insert(con, "hsp", handler$getHsp()))
  new_blastReportDB(con)
}


#' @usage blastReportDBConnect(db_path)
#' @return A \code{\linkS4class{blastReportDB}} object.
#' @rdname blastReportDB
#' @export
blastReportDBConnect <- function(db_path) {
  assert_that(is.readable(db_path))
  new_blastReportDB(db_connect(db_path))
}
