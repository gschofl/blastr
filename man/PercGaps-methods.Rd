\name{getPercGaps}
\alias{getPercGaps}
\alias{getPercGaps,Hit-method}
\alias{getPercGaps,HitList-method}
\alias{getPercGaps,Hsp-method}
\alias{getPercGaps,HspList-method}
\alias{getPercGaps,Iteration-method}
\alias{getPercGaps-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getPercGaps(x, max = FALSE)

  \S4method{getPercGaps}{Hsp}(x, ...)

  \S4method{getPercGaps}{HspList}(x, ...)

  \S4method{getPercGaps}{Hit}(x, ...)

  \S4method{getPercGaps}{HitList}(x, ...)

  \S4method{getPercGaps}{Iteration}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
  An numeric vector or a list of numeric vectors
}
\description{
  Extract Percent Gaps.
}

