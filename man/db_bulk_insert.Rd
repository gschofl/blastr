\docType{methods}
\name{db_bulk_insert}
\alias{db_bulk_insert}
\title{Insert a \code{data.frame} into a corresponding SQLite table}
\arguments{
  \item{con}{A connection object.}

  \item{tbl}{Name of table in database.}

  \item{df}{A \code{data.frame} matching \code{tbl}.}

  \item{...}{Additional arguments (Currently \code{log}:
  Path to a log file).}
}
\description{
  Insert a \code{data.frame} into a corresponding SQLite
  table
}
\keyword{internal}

