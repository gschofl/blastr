\name{db_create}
\alias{db_create}
\title{Create an SQLite database.}
\usage{
  db_create(dbName, dbSchema = "", overwrite = FALSE)
}
\arguments{
  \item{dbName}{Path to an SQLite database.}

  \item{dbSchema}{SQL schema for setting up the db.}

  \item{overwrtite}{Overwrite an existing db file by the
  same name.}
}
\description{
  Create an SQLite database.
}
\keyword{internal}

