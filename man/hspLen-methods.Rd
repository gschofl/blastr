\name{hspLen}
\alias{hspLen}
\alias{hspLen,Hit-method}
\alias{hspLen,HitList-method}
\alias{hspLen,Iteration-method}
\alias{hspLen,IterationList-method}
\alias{hspLen-methods}
\title{Access components of Blast Hits}
\usage{
  hspLen(x)

  \S4method{hspLen}{Hit}(x, ...)

  \S4method{hspLen}{HitList}(x, ...)

  \S4method{hspLen}{Iteration}(x, ...)

  \S4method{hspLen}{IterationList}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}
}
\value{
  An numeric vector.
}
\description{
  Extract the number of Hsps for each Hit.
}

