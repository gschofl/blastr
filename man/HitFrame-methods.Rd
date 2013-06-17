\name{getHitFrame}
\alias{getHitFrame}
\alias{getHitFrame,Hit-method}
\alias{getHitFrame,HitList-method}
\alias{getHitFrame,Hsp-method}
\alias{getHitFrame,HspList-method}
\alias{getHitFrame,Iteration-method}
\alias{getHitFrame,blastReportDB-method}
\alias{getHitFrame-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getHitFrame(x, max = FALSE)

  \S4method{getHitFrame}{Hsp}(x, ...)

  \S4method{getHitFrame}{HspList}(x, ...)

  \S4method{getHitFrame}{Hit}(x, ...)

  \S4method{getHitFrame}{HitList}(x, ...)

  \S4method{getHitFrame}{Iteration}(x, ...)

  \S4method{getHitFrame}{blastReportDB}(x, ...)
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
  Extract the frame of the hit sequence.
}

