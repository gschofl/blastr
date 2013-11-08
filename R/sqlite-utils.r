#' @include utils.r
#' @importClassesFrom RSQLite SQLiteConnection
#' @importClassesFrom RSQLite SQLiteObject
#' @importClassesFrom RSQLite dbObjectId
#' @importClassesFrom DBI DBIConnection
#' @importClassesFrom DBI DBIObject
#' @importFrom RSQLite dbListTables
#' @importFrom RSQLite SQLite
#' @importFrom RSQLite dbConnect
#' @importFrom RSQLite dbGetQuery
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite dbCommit
#' @importFrom RSQLite dbListTables
#' @importFrom RSQLite dbListFields
#' @importFrom RSQLite dbBeginTransaction
#' @importFrom RSQLite dbSendPreparedQuery
#' @importFrom RSQLite dbGetInfo
NULL

#' Create an SQLite database.
#' 
#' @param dbName Path to an SQLite database.
#' @param dbSchema SQL schema for setting up the db.
#' @keywords internal
db_create <- function(dbName, dbSchema = "") {
  assert_that(is.string(dbSchema))
  message('Creating new database ', sQuote(basename(dbName)))
  if (file.exists(dbName)) {
    unlink(dbName)
  }
  con <- dbConnect(SQLite(), dbname = dbName)
  sql <- compactChar(trim(strsplit(dbSchema, ";\n")[[1L]]))
  if (length(sql) > 0L) {
    tryCatch(lapply(sql, dbGetQuery, conn = con), error = function(e) {
      message(e)
    })
  } 
  con
}


#' Connect to an existing SQLite database.
#' 
#' @param dbName Path to an SQLite database.
#' @param message Message if db does not exits
#' @keywords internal
db_connect <- function(dbName, message = "") {
  if (!file.exists(dbName))
    stop("Database ", sQuote(basename(dbName)),
         " does not exist.\n", message, call.=FALSE)
  dbConnect(SQLite(), dbname=dbName)
}


#' Disconnect from an SQLite database.
#' 
#' @param ... connection objects.
#' @keywords internal
db_disconnect <- function(...) {
  lapply(list(...), dbDisconnect)
}


#' Query an SQLite database.
#' 
#' @param con a connection object.
#' @param sql an SQL statemant
#' @param j
#' @keywords internal
db_query <- function(con, sql, j=NA) {
  assert_that(is(con, "SQLiteConnection"))
  assert_that(is.string(sql), noNA(sql))
  data <- dbGetQuery(con, sql)
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
#' @keywords internal
db_count <- function(con, tbl) {
  assert_that(con %has_tables% tbl)
  sql <- paste0("SELECT count(*) FROM ", tbl)
  db_query(con, sql, 1)
}


#' @keywords internal
db_bulk_insert <- function(con, tbl, df) {
  sql <- sprintf("INSERT INTO %s VALUES (%s)", tbl,
                 paste0("$", names(df), collapse=", "))
  dbBeginTransaction(con)
  dbSendPreparedQuery(con, sql, df)
  dbCommit(con)
}


#' Check if a SQLite database has specified tables
#'
#' @param con a connection object
#' @param tables a character vector of table names
#' @keywords internal
has_tables <- function (con, tbl) {
  assert_that(is(con, "SQLiteConnection"))
  all(tbl %in% dbListTables(con))
}
on_failure(has_tables) <- function(call, env) {
  tbl <- paste0(eval(call$tbl, env), collapse = ", ")
  dbName <- dbGetInfo(eval(call$con, env))$dbname
  paste0("Missing table(s) ", tbl, " in database ", sQuote(dbName))
}


#' @keywords internal
#' @rdname has_tables
"%has_tables%" <- has_tables

