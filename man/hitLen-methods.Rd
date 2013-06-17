\name{hitLen}
\alias{hitLen}
\alias{hitLen,Iteration-method}
\alias{hitLen,IterationList-method}
\alias{hitLen-methods}
\alias{hspLen,IterationList-method}
\title{Access components of Blast Iterations (Queries)}
\usage{
  hitLen(x)

  \S4method{hitLen}{Iteration}(x, ...)

  \S4method{hitLen}{IterationList}(x, ...)
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

