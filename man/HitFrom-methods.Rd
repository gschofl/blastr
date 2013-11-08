\name{getHitFrom}
\alias{getHitFrom}
\alias{getHitFrom,Hit-method}
\alias{getHitFrom,HitList-method}
\alias{getHitFrom,Hsp-method}
\alias{getHitFrom,HspList-method}
\alias{getHitFrom,Iteration-method}
\alias{getHitFrom,blastReportDB-method}
\alias{getHitFrom-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getHitFrom(x, max = FALSE)

  \S4method{getHitFrom}{Hsp}(x, ...)

  \S4method{getHitFrom}{HspList}(x, ...)

  \S4method{getHitFrom}{Hit}(x, ...)

  \S4method{getHitFrom}{HitList}(x, ...)

  \S4method{getHitFrom}{Iteration}(x, ...)

  \S4method{getHitFrom}{blastReportDB}(x, ...)
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
  Extract the position of the start residues for the hit
  sequence.
}

