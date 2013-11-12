\name{db_connect}
\alias{db_connect}
\title{Connect to an existing SQLite database.}
\usage{
  db_connect(dbName, message = "", create = FALSE)
}
\arguments{
  \item{dbName}{Path to an SQLite database.}

  \item{message}{Message to be produced if db does not
  exist.}

  \item{create}{If \code{TRUE}, create the database if it
  doesn't exist.}
}
\description{
  Connect to an existing SQLite database.
}
\keyword{internal}

