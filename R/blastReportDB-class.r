#' @include blastReport-class.r
#' @include blast-parser-db.r
#' @include sqlite-db.r
#' @importFrom IRanges IRanges IRangesList reduce width unlist
NULL


# blastReportDB-class ----------------------------------------------------


#' blastReportDB-class
#' 
#' \sQuote{\bold{blastReportDB}}: A connection to an SQLite database
#' containing blast records organised in three tables:
#' 
#' \bold{query} with fields:
#' 
#' \itemize{
#'    \item query_id     INTEGER    Primary key
#'    \item query_def    TEXT
#'    \item query_len    INTEGER
#' }
#' 
#' \bold{hit} with fields:
#' 
#' \itemize{
#'    \item query_id     INTEGER
#'    \item hit_id       INTEGER    Primary key
#'    \item hit_num      INTEGER
#'    \item gene_id      TEXT
#'    \item accession    TEXT
#'    \item definition   TEXT
#'    \item length       INTEGER
#' }
#' 
#' \bold{hsp} with fields:
#'
#' \itemize{
#'    \item query_id     INTEGER
#'    \item hit_id       INTEGER    
#'    \item hsp_id       INTEGER    Primary key
#'    \item hsp_num      INTEGER
#'    \item bit_score    FLOAT
#'    \item score        INTEGER
#'    \item evalue       FLOAT
#'    \item query_from   INTEGER
#'    \item query_to     INTEGER
#'    \item hit_from     INTEGER
#'    \item hit_to       INTEGER
#'    \item query_frame  INTEGER
#'    \item query_frame  INTEGER
#'    \item identity     INTEGER
#'    \item positive     INTEGER
#'    \item gaps         INTEGER
#'    \item align_len    INTEGER
#'    \item qseq         TEXT
#'    \item hseq         TEXT
#'    \item midline      TEXT
#' }
#'    
#' @seealso
#'  The constructors \code{\link{blastReportDB}} and \code{\link{blastReportDBConnect}},
#'  and the BLAST classes \code{\linkS4class{blastReport}} and
#'  \code{\linkS4class{blastTable}}.
#'   
#' @name blastReportDB-class
#' @rdname blastReportDB-class
#' @exportClass blastReportDB
.blastReportDB <- setRefClass(
  Class='blastReportDB',
  contains='sqliteDB',
  methods=list(
    initialize=function(con = db_create(dbSchema=blast_db.sql(), verbose=FALSE), ...) {
      callSuper(con, ...)
    })
)


setValidity('blastReportDB', function(object) {
  errors <- character()
  if (length(db_list_tables(object)) == 0L) {
    return("No tables in 'blastReportDB'")
  }
  if (!all(c("hit", "hsp", "query") %in% db_list_tables(object))) {
    errors <- c(errors, "Table missing from 'blastReportDB'\n")
  }
  if (!all(c("query_id", "query_def", "query_len") %in% db_list_fields(object, "query"))) {
    errors <- c(errors, "Field missing from table 'query'\n")
  }
  if (!all(c("query_id", "hit_id", "hit_num", "gene_id", "accession",
             "definition", "length") %in% db_list_fields(object, "hit"))) {
    errors <- c(errors, "Field missing from table 'hit'\n")
  }
  if (!all(c("query_id", "hit_id", "hsp_id", "hsp_num", "bit_score",
             "score", "evalue", "query_from", "query_to", "hit_from",
             "hit_to", "query_frame", "hit_frame", "identity", "positive",
             "gaps", "align_len", "qseq", "hseq", "midline") 
           %in% db_list_fields(object, "hsp"))) {
    errors <- c(errors, "Field missing from table 'hsp'")
  }
  
  if (length(errors) == 0L) TRUE else errors
})


#' @aliases show,blastReportDB-method
#' @rdname show-methods
setMethod('show', 'blastReportDB',
          function(object) {
            showme <- sprintf('%s connection object with:\n| %s queries | %s hits | %s hsps |',
                              sQuote(class(object)),
                              db_count(object, "query"),
                              db_count(object, "hit"),
                              db_count(object, "hsp"))
            cat(showme, sep="\n")
          })


setMethod('lapply', 'blastReportDB', function(X, FUN, ...) {
  lapply(db_query(X, "select query_id from query", 1L), function(j)
    FUN(X[j], ...)
  )
})


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
#' @param verbose Message if a new database is created.
#'
#' @return A \code{\linkS4class{blastReportDB}} object.
#' @rdname blastReportDB
#' @export
blastReportDB <- function(blastfile, db_path = ":memory:", max_hit = NULL,
                          max_hsp = NULL, reset_at = 1000, verbose = TRUE) {
  db <- .blastReportDB(db_create(db_path, blast_db.sql(), overwrite=TRUE, verbose=verbose))
  if (missing(blastfile)) {
    return(db)
  }
  assert_that(is.readable(blastfile), has_extension(blastfile, 'xml'))
  handler <- BlastOutput.Iterations(conn(db), max_hit, max_hsp, reset_at)
  out <- xmlEventParse(blastfile, list(), branches=handler)
  ## load final part into db
  assert_that(.db_bulk_insert(con=conn(db), tbl="query", df=handler$getQuery()))
  assert_that(.db_bulk_insert(conn(db), "hit", handler$getHit()))
  assert_that(.db_bulk_insert(conn(db), "hsp", handler$getHsp()))
  validObject(db)
  db
}


#' @usage blastReportDBConnect(db_path)
#' @return A \code{\linkS4class{blastReportDB}} object.
#' @rdname blastReportDB
#' @export
blastReportDBConnect <- function(db_path) {
  if (db_path == ":memory:") {
    stop("Cannot connect to an in-memory database", call.=FALSE)
  }
  if (db_path == "") {
    stop("Cannot connect to a temporary database", call.=FALSE)
  }
  db <- .blastReportDB(db_connect(db_path))
  validObject(db)
  db
}


