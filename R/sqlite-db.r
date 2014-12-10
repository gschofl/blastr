#' @include utils.r
#' @importFrom assertthat is.readable
#' @importClassesFrom RSQLite SQLiteConnection
#' @importClassesFrom DBI DBIConnection DBIObject
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbCommit 
#' @importFrom RSQLite dbGetQuery dbBegin dbSendPreparedQuery
#' @importFrom RSQLite dbListTables dbListFields dbGetInfo dbIsValid
NULL


# SQLite utils ------------------------------------------------------------


#' Connect to an existing SQLite database.
#' 
#' @param dbName Path to an SQLite database [default: :memory:].
#' @param create
#' @keywords internal
db_connect <- function(dbName = ":memory:", create = FALSE) {
  assert_that(is.string(dbName))
  if (dbName != ":memory:" && dbName != "" && !create) {
    assert_that(is.readable(dbName))
  }
  dbConnect(SQLite(), dbname = dbName)
}

#' Disconnect from an SQLite database.
#' 
#' @param ... Connection objects.
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
#' @keywords internal
db_create <- function(dbName = ":memory:", dbSchema = "", overwrite = FALSE, verbose = TRUE) {
  assert_that(is.string(dbName), is.string(dbSchema))
  if (file.exists(dbName)) {
    if (overwrite) {
      unlink(dbName)
    } else {
      stop("File ", sQuote(basename(dbName)), " already exists. Use 'db_connect'.", call. = FALSE)  
    }  
  }
  if (verbose) {
    message('Creating new database ', sQuote(basename(dbName)))
  }
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
#' @keywords internal
db_info <- function(con) {
  dbGetInfo(con)
}

#' Is a connection to an SQLite database still current.
#' 
#' @param con A connection object.
#' @keywords internal
db_is_current <- function(con) {
  dbIsValid(con)
}

#' Check if a SQLite database has specified tables
#'
#' @param con a connection object
#' @param tables a character vector of table names
#' @param ...
#' @keywords internal
has_tables <- function(con, tbl, ...) {
  all(tbl %in% db_list_tables(con, ...))
}
on_failure(has_tables) <- function(call, env) {
  tbl <- paste0(eval(call$tbl, env), collapse = ", ")
  dbName <- dbGetInfo(eval(call$con, env))$dbname
  paste0("Missing table(s) ", tbl, " in database ", sQuote(dbName))
}

#' @keywords internal
"%has_tables%" <- has_tables


# SQLiteDB-class ---------------------------------------------------------


#' Class \code{"SQLiteDB"}
#' 
#' @field .con An \code{\linkS4class{SQLiteConnection}}.
#' @field .info A list.
#' @field .path Path to the database.
#' @keywords internal   
sqliteDB <- setRefClass(
  Class = 'SQLiteDB',
  fields = list(
    .con = 'SQLiteConnection',
    .info = 'list',
    .path = 'character'
  ),
  methods = list(
    initialize = function(con = dbConnect(RSQLite::SQLite(), dbname = ":memory:"), ...) {
      .con <<- con
      .info <<- dbGetInfo(con)
      .path <<- if (is.null(.info$dbname)) {
          ''
        } else {
          normalizePath(.info$dbname, mustWork = TRUE)
        }
    })
)

#' @keywords internal
setGeneric("conn", function(x) standardGeneric("conn"))
setMethod("conn", "SQLiteDB", function(x) x$.con)

#' @keywords internal
setGeneric("info", function(x) standardGeneric("info"))
setMethod("info", "SQLiteDB", function(x) x$.info)

#' @keywords internal
setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "SQLiteDB", function(x) x$.path)


# SQLite transactions -----------------------------------------------------


#' Query an SQLite database.
#' 
#' @usage db_query(x, stmt, j = NA, ...)
#' @param x A connection object.
#' @param stmt An SQL statemant.
#' @param j Subset the result.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @keywords internal
setGeneric("db_query", function(x, stmt, j = NA, ...) standardGeneric("db_query"))

.db_query <- function(con, stmt, j = NA, log = NULL) {  
  do_log(ellipsize(stmt, 60), '\n', to = log)
  data <- tryCatch(dbGetQuery(con, stmt), error = function(e) {
    do_log(e$message, '\n', to = log)
    return(NULL)
  })
  if (is.na(j)) {
    return(data)
  } else if (nrow(data) == 0) {
    return(character(0))
  } else {
    return(data[[j]])
  }
}

setMethod("db_query", "SQLiteDB", function(x, stmt, j = NA, ...) {
  .db_query(conn(x), stmt = stmt, j = j, log = list(...)$log)
})

setMethod("db_query", "SQLiteConnection", function(x, stmt, j = NA, ...) {
  .db_query(x, stmt = stmt, j = j, log = list(...)$log)
})

#' Insert a \code{data.frame} into a corresponding SQLite table
#' 
#' @usage db_bulk_insert(x, tbl, df, ...)
#' @param con A connection object.
#' @param tbl Name of table in database.
#' @param df A \code{data.frame} matching \code{tbl}.
#' @param ... Additional arguments (Currently \code{log}: Path to a log file). 
#' @keywords internal
setGeneric("db_bulk_insert", function(x, tbl, df, ...) standardGeneric("db_bulk_insert"))

.db_bulk_insert <- function(con, tbl, df, log = NULL) {
  stmt <- sprintf("insert into %s values (%s)", tbl, comma("$", names(df)))
  dbBegin(con)
  tryCatch(dbSendPreparedQuery(con, stmt, df), error = function(e) {
    do_log(e$message, '\n', to = log)
    return(NULL)
  })
  dbCommit(con)
}

setMethod("db_bulk_insert", "SQLiteDB", function(x, tbl, df, ...) {
  .db_bulk_insert(conn(x), tbl, df, log = list(...)$log)
})

setMethod("db_bulk_insert", "SQLiteConnection", function(x, tbl, df, ...) {
  .db_bulk_insert(x, tbl, df, log = list(...)$log)
})

#' List available tables in an SQLite database
#' 
#' @usage db_list_tables(x, ...)
#' @param x A connection object.
#' @param ... Additional arguments (Currently unused). 
#' @keywords internal
setGeneric("db_list_tables", function(x, ...) standardGeneric("db_list_tables"))

setMethod("db_list_tables", "SQLiteDB", function(x, ...) dbListTables(conn(x)))

setMethod("db_list_tables", "SQLiteConnection", function(x, ...) dbListTables(x))

#' List available tables in an SQLite database
#' 
#' @usage db_list_fields(x, tbl, ...)
#' @param x A connection object.
#' @param tbl Name of the table.
#' @param ... Additional arguments (Currently unused). 
#' @keywords internal
setGeneric("db_list_fields", function(x, tbl, ...) standardGeneric("db_list_fields"))

setMethod("db_list_fields", "SQLiteDB", function(x, tbl, ...) dbListFields(conn(x), tbl))

setMethod("db_list_fields", "SQLiteConnection", function(x, tbl, ...) dbListFields(x, tbl))

#' Count rows in a db table
#' 
#' @usage db_count(x, tbl)
#' @param x A connection object.
#' @param tbl Name of table in database.
#' @param ... Additional arguments (Currently unused). 
#' @keywords internal
setGeneric("db_count", function(x, tbl, ...) standardGeneric("db_count"))

.db_count <- function(con, tbl, ...) {
  assert_that(has_tables(con, tbl, ...))
  .db_query(con, paste0("select count(*) from ", tbl), 1L, ...)
}

setMethod("db_count", "SQLiteDB", function(x, tbl, ...) .db_count(conn(x), tbl, ...))

setMethod("db_count", "SQLiteConnection", function(x, tbl, ...) .db_count(x, tbl, ...))

