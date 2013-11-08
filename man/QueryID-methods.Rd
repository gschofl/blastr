\docType{methods}
\name{getQueryID}
\alias{getQueryID}
\alias{getQueryID,Iteration-method}
\alias{getQueryID,IterationList-method}
\alias{getQueryID,blastReport-method}
\alias{getQueryID,blastReportDB-method}
\alias{getQueryID-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  getQueryID(x, ...)

  \S4method{getQueryID}{Iteration}(x, ...)

  \S4method{getQueryID}{IterationList}(x, ...)

  \S4method{getQueryID}{blastReport}(x, ...)

  \S4method{getQueryID}{blastReportDB}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{blastReport} or
  \linkS4class{IterationList} object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  A numeric vector or list of numeric vectors.
}
\description{
  Extract the query ID.
}

