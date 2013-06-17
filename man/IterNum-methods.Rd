\name{getIterNum}
\alias{getIterNum}
\alias{getIterNum,Iteration-method}
\alias{getIterNum,IterationList-method}
\alias{getIterNum,blastReport-method}
\alias{getIterNum-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  getIterNum(x, ...)

  \S4method{getIterNum}{Iteration}(x, ...)

  \S4method{getIterNum}{IterationList}(x, ...)

  \S4method{getIterNum}{blastReport}(x, ...)
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
  Extract the Iteration (Query) numbers.
}

