#' @include utils.r
#' @importFrom assertthat is.readable
#' @importClassesFrom RSQLite SQLiteConnection SQLiteObject dbObjectId
#' @importClassesFrom DBI DBIConnection DBIObject
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbCommit 
#' @importFrom RSQLite dbGetQuery dbBeginTransaction sqliteExecStatement
#' @importFrom RSQLite dbListTables dbListFields dbGetInfo isIdCurrent
NULL

#' sqliteDB-class
#' 
#' @keywords internal   
#' @name sqliteDB-class
#' @rdname sqliteDB-class
#' @exportClass sqliteDB
.sqliteDB <- setRefClass(
  Class='sqliteDB',
  fields=list(
    .con = 'SQLiteConnection',
    .info = 'list',
    .path = 'character'
  ),
  methods=list(
    initialize=function(con = db_create(':memory:', verbose=FALSE), ...) {
      .con <<- con
      .info <<- db_info(.con)
      .path <<- 
        if (.info$dbname == ':memory:' || .info$dbname == '') {
          ''
        } else {
          normalizePath(.info$dbname, mustWork=TRUE)
        }
    })
)

now <- function(accuracy = 4) {
  paste0("-- ", format(Sys.time(), paste0("%M:%OS", accuracy)), " -- ")
}

do_log <- function(log = NULL, ...) {
  if (!is.null(log)) {
    cat(now(), ..., file=log, append=TRUE, sep="")
  }
}


#' @keywords internal
#' @export
#' @docType methods
setGeneric("conn", function(x) standardGeneric("conn"))
setMethod("conn", "sqliteDB", function(x) x$.con)

#' @keywords internal
#' @export
#' @docType methods
setGeneric("info", function(x) standardGeneric("info"))
setMethod("info", "sqliteDB", function(x) x$.info)

#' @keywords internal
#' @export
#' @docType methods
setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "sqliteDB", function(x) x$.path)


.db_query <- function(con, stmt, j = NA, log = NULL) {
  do_log(log, ellipsize(stmt, 60), '\n')
  data <- tryCatch(dbGetQuery(con, stmt), error = function(e) {
    do_log(log, e$message, '\n')
    return(NULL)
  })
  if (is.na(j))
    return(data)
  else if (nrow(data) == 0)
    return(character(0))
  else
    return(data[[j]])
}
#' Query an SQLite database.
#' 
#' @usage db_query(x, stmt, j = NA, ...)
#' @param x A connection object.
#' @param stmt An SQL statemant.
#' @param j Subset the result.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("db_query", function(x, stmt, j = NA, ...) standardGeneric("db_query"))
setMethod("db_query", "sqliteDB", function(x, stmt, j = NA, ...) {
  log <- list(...)$log
  .db_query(conn(x), stmt = stmt, j = j, log = log)
})
setMethod("db_query", "SQLiteConnection", function(x, stmt, j = NA, ...) {
  log <- list(...)$log
  .db_query(x, stmt = stmt, j = j, log = log)
})


.db_bulk_insert <- function(con, tbl, df, log = NULL) {
  stmt <- sprintf("insert into %s values (%s)", tbl, paste0("$", names(df), collapse=", "))
  dbBeginTransaction(con)
  do_log(log, 'Exec statement [', tbl, '] ')
  tryCatch(sqliteExecStatement(con, stmt, df), error = function(e) {
    do_log(log, e$message, '\n')
    return(NULL)
  })
  rc <- dbCommit(con)
  do_log(log, 'Commit [', tbl, '] ', rc, '\n')
  rc
}
#' Insert a \code{data.frame} into a corresponding SQLite table
#' 
#' @usage db_bulk_insert(x, tbl, df, ...)
#' @param con A connection object.
#' @param tbl Name of table in database.
#' @param df A \code{data.frame} matching \code{tbl}.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("db_bulk_insert", function(x, tbl, df, ...) standardGeneric("db_bulk_insert"))
setMethod("db_bulk_insert", "sqliteDB", function(x, tbl, df, ...) {
  log <- list(...)$log
  .db_bulk_insert(conn(x), tbl, df, log)
})
setMethod("db_bulk_insert", "SQLiteConnection", function(x, tbl, df, ...) {
  log <- list(...)$log
  .db_bulk_insert(x, tbl, df, log)
})


#' List available tables in an SQLite database
#' 
#' @usage db_list_tables(x, ...)
#' @param x A connection object.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("db_list_tables", function(x, ...) standardGeneric("db_list_tables"))
setMethod("db_list_tables", "sqliteDB", function(x, ...) {
  log <- list(...)$log
  do_log(log, 'list tables\n')
  dbListTables(conn(x))
})
setMethod("db_list_tables", "SQLiteConnection", function(x, ...) {
  log <- list(...)$log
  do_log(log, 'list tables\n')
  dbListTables(x)
})


#' List available tables in an SQLite database
#' 
#' @usage db_list_fields(x, tbl, ...)
#' @param x A connection object.
#' @param tbl Name of the table.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("db_list_fields", function(x, tbl, ...) standardGeneric("db_list_fields"))
setMethod("db_list_fields", "sqliteDB", function(x, tbl, ...) {
  log <- list(...)$log
  do_log(log, 'list fields in ', tbl, '\n')
  dbListFields(conn(x), tbl)
})
setMethod("db_list_fields", "SQLiteConnection", function(x, tbl, ...) {
  log <- list(...)$log
  do_log(log, 'list fields in ', tbl, '\n')
  dbListFields(x, tbl)
})


.db_count <- function(con, tbl, ...) {
  assert_that(has_tables(con, tbl, ...))
  stmt <- paste0("select count(*) from ", tbl)
  db_query(con, stmt, 1L, ...)
}
#' Count rows in a db table
#' 
#' @usage db_count(x, tbl)
#' @param x A connection object.
#' @param tbl Name of table in database.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("db_count", function(x, tbl, ...) standardGeneric("db_count"))
setMethod("db_count", "sqliteDB", function(x, tbl, ...) {
  .db_count(x, tbl, ...)
})
setMethod("db_count", "SQLiteConnection", function(x, tbl, ...) {
  .db_count(x, tbl, ...)
})


# SQLite utils ------------------------------------------------------------



#' Connect to an existing SQLite database.
#' 
#' @param dbName Path to an SQLite database [default: :memory:].
#' @param create
#' @export
#' @keywords internal
db_connect <- function(dbName = ":memory:", create = FALSE) {
  assert_that(is.string(dbName))
  if (dbName != ":memory:" && dbName != "" && !create) {
    assert_that(is.readable(dbName))
  }
  con <- dbConnect(SQLite(), dbname = dbName)
  con
}


#' Disconnect from an SQLite database.
#' 
#' @param ... Connection objects.
#' @export
#' @keywords internal
db_disconnect <- function(...) {
  lapply(list(...), dbDisconnect)
}


#' Create an SQLite database.
#' 
#' @param dbName Path to an SQLite database.
#' @param dbSchema SQL schema for setting up the db.
#' @param overwrite Overwrite an existing db file by the same name.
#' @param verbose
#' @export
#' @keywords internal
db_create <- function(dbName = ":memory:", dbSchema = "", overwrite = FALSE, verbose = TRUE) {
  assert_that(is.string(dbName), is.string(dbSchema))
  if (file.exists(dbName)) {
    if (overwrite)
      unlink(dbName)
    else
      stop("File ", sQuote(basename(dbName)), " already exists. Use 'db_connect'.", call.=FALSE)
  }
  if (verbose)
    message('Creating new database ', sQuote(basename(dbName)))
  con <- db_connect(dbName = dbName, create = TRUE)
  sql <- compactChar(trim(strsplit(dbSchema, ";\n")[[1L]]))
  if (length(sql) > 0L) {
    tryCatch(lapply(sql, dbGetQuery, conn = con), error = function(e) {
      message(e)
    })
  } 
  con
}


#' Metadata for an SQLite database.
#' 
#' @param con A connection object.
#' @export
#' @keywords internal
db_info <- function(con) {
  dbGetInfo(con)
}


#' Is a connection to an SQLite database still current.
#' 
#' @param con A connection object.
#' @export
#' @keywords internal
db_is_current <- function(con) {
  isIdCurrent(con)
}


#' Check if a SQLite database has specified tables
#'
#' @param con a connection object
#' @param tables a character vector of table names
#' @param ...
#' @export
#' @keywords internal
has_tables <- function(con, tbl, ...) {
  all(tbl %in% db_list_tables(con, ...))
}
on_failure(has_tables) <- function(call, env) {
  tbl <- paste0(eval(call$tbl, env), collapse = ", ")
  dbName <- dbGetInfo(eval(call$con, env))$dbname
  paste0("Missing table(s) ", tbl, " in database ", sQuote(dbName))
}

#' @export
#' @keywords internal
"%has_tables%" <- has_tables
