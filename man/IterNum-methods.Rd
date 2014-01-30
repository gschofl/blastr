\docType{methods}
\name{getIterNum,Iteration-method}
\alias{getIterNum}
\alias{getIterNum,Iteration-method}
\alias{getIterNum,IterationList-method}
\alias{getIterNum,blastReport-method}
\title{Access components of Blast Iterations (Queries)}
\usage{
\S4method{getIterNum}{Iteration}(x)

\S4method{getIterNum}{IterationList}(x)

getIterNum(x, ...)

\S4method{getIterNum}{blastReport}(x)
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

