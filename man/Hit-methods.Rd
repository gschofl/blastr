\docType{methods}
\name{getHit,Iteration-method}
\alias{getHit}
\alias{getHit,Iteration-method}
\alias{getHit,IterationList-method}
\alias{getHit,blastReport-method}
\title{Access components of Blast Iterations (Queries)}
\usage{
\S4method{getHit}{Iteration}(x, n = NULL, drop = TRUE, ...)

\S4method{getHit}{IterationList}(x, n = NULL, drop = TRUE, ...)

getHit(x, ...)

\S4method{getHit}{blastReport}(x, n = NULL, drop = TRUE)
}
\arguments{
  \item{x}{A \linkS4class{blastReport} or
  \linkS4class{IterationList} object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
A \linkS4class{hit} or \linkS4class{HitList} object.
}
\description{
Extract Blast hits.
}

