\name{db_create}
\alias{db_create}
\title{Create an SQLite database.}
\usage{
  db_create(dbName = ":memory:", dbSchema = "",
    overwrite = FALSE, verbose = TRUE)
}
\arguments{
  \item{dbName}{Path to an SQLite database.}

  \item{dbSchema}{SQL schema for setting up the db.}

  \item{overwrite}{Overwrite an existing db file by the
  same name.}

  \item{verbose}{}
}
\description{
  Create an SQLite database.
}
\keyword{internal}

