#' @include utils.r
#' @importFrom assertthat assert_that is.string is.readable not_empty noNA
#' @importClassesFrom RSQLite SQLiteConnection SQLiteObject dbObjectId
#' @importClassesFrom DBI DBIConnection DBIObject
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbCommit 
#' @importFrom RSQLite dbGetQuery dbBeginTransaction dbSendPreparedQuery
#' @importFrom RSQLite dbListTables dbListFields dbGetInfo
NULL


#' Connect to an existing SQLite database.
#' 
#' @param dbName Path to an SQLite database.
#' @param message Message to be produced if db does not exist.
#' @param create If \code{TRUE}, create the database if it doesn't exist.
#' @export
#' @keywords internal
db_connect <- function(dbName, message = "", create = FALSE) {
  if (!file.exists(dbName) && !create)
    stop("Database ", sQuote(basename(dbName)), " does not exist.\n", message, call.=FALSE)
  dbConnect(SQLite(), dbname=dbName)
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
#' @param overwrtite Overwrite an existing db file by the same name.
#' @export
#' @keywords internal
db_create <- function(dbName, dbSchema = "", overwrite = FALSE) {
  assert_that(is.string(dbName), is.string(dbSchema))
  if (file.exists(dbName)) {
    if (overwrite)
      unlink(dbName)
    else
      stop("File ", sQuote(basename(dbName)), " already exists. Use 'db_connect'.", call.=FALSE)
  }
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

#' Query an SQLite database.
#' 
#' @param con a connection object.
#' @param stmt an SQL statemant.
#' @param j Subset data.
#' @export
#' @keywords internal
db_query <- function(con, stmt, j=NA) {
  assert_that(is(con, "SQLiteConnection"), is.string(stmt), noNA(stmt))
  data <- dbGetQuery(con, stmt)
  if (is.na(j))
    return(data)
  if (nrow(data) == 0)
    return(character(0))
  else
    return(data[[j]])
}

#' Count rows in a db table
#' 
#' @param con a connection object.
#' @param tbl name of table in database.
#' @export
#' @keywords internal
db_count <- function(con, tbl) {
  assert_that(con %has_tables% tbl)
  sql <- paste0("select count(*) from ", tbl)
  db_query(con, sql, 1)
}


#' Bulk insert a data.frame into a db table
#' 
#' @param con A connection object.
#' @param tbl Name of table in database.
#' @param df A \code{data.frame} matching the \code{tbl}.
#' @export
#' @keywords internal
db_bulk_insert <- function(con, tbl, df) {
  sql <- sprintf("insert into %s values (%s)", tbl,
                 paste0("$", names(df), collapse=", "))
  dbBeginTransaction(con)
  dbSendPreparedQuery(con, sql, df)
  dbCommit(con)
}

#' Check if a SQLite database has specified tables
#'
#' @param con a connection object
#' @param tables a character vector of table names
#' @export
#' @keywords internal
has_tables <- function(con, tbl) {
  assert_that(is(con, "SQLiteConnection"))
  all(tbl %in% dbListTables(con))
}
on_failure(has_tables) <- function(call, env) {
  tbl <- paste0(eval(call$tbl, env), collapse = ", ")
  dbName <- dbGetInfo(eval(call$con, env))$dbname
  paste0("Missing table(s) ", tbl, " in database ", sQuote(dbName))
}

#' @export
#' @keywords internal
"%has_tables%" <- has_tables

