\name{getHit}
\alias{getHit}
\alias{getHit,Iteration-method}
\alias{getHit,IterationList-method}
\alias{getHit,blastReport-method}
\alias{getHit-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  getHit(x, ...)

  \S4method{getHit}{Iteration}(x, ...)

  \S4method{getHit}{IterationList}(x, ...)

  \S4method{getHit}{blastReport}(x, ...)
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

