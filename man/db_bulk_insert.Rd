\name{db_bulk_insert}
\alias{db_bulk_insert}
\title{Bulk insert a data.frame into a db table}
\usage{
  db_bulk_insert(con, tbl, df)
}
\arguments{
  \item{con}{A connection object.}

  \item{tbl}{Name of table in database.}

  \item{df}{A \code{data.frame} matching the \code{tbl}.}
}
\description{
  Bulk insert a data.frame into a db table
}
\keyword{internal}

