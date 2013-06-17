\name{getIteration}
\alias{getIteration}
\alias{getIteration,blastReport-method}
\alias{getIteration-methods}
\title{Extract Iterations/Queries from a BLAST report}
\usage{
  getIteration(x, ...)

  \S4method{getIteration}{blastReport}(x, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{blastReport}} object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  A \linkS4class{Iteration} or \linkS4class{IterationList}
  object.
}
\description{
  Extract Iterations/Queries from a BLAST report
}

