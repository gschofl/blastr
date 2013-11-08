\name{nhits}
\alias{nhits}
\alias{nhits,Iteration-method}
\alias{nhits,IterationList-method}
\alias{nhits-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  nhits(x)

  \S4method{nhits}{Iteration}(x, ...)

  \S4method{nhits}{IterationList}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{blastReport} or
  \linkS4class{IterationList} object.}
}
\value{
  An numeric vector.
}
\description{
  Extract the number of Hits for each Iteration/Query.
}

