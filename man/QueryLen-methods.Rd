\name{getQueryLen}
\alias{getQueryLen}
\alias{getQueryLen,Iteration-method}
\alias{getQueryLen,IterationList-method}
\alias{getQueryLen,blastReport-method}
\alias{getQueryLen,blastReportDB-method}
\alias{getQueryLen-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  getQueryLen(x, ...)

  \S4method{getQueryLen}{Iteration}(x, ...)

  \S4method{getQueryLen}{IterationList}(x, ...)

  \S4method{getQueryLen}{blastReport}(x, ...)

  \S4method{getQueryLen}{blastReportDB}(x, ...)
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
  Extract the query length.
}

