\name{getQueryFrame}
\alias{getQueryFrame}
\alias{getQueryFrame,Hit-method}
\alias{getQueryFrame,HitList-method}
\alias{getQueryFrame,Hsp-method}
\alias{getQueryFrame,HspList-method}
\alias{getQueryFrame,Iteration-method}
\alias{getQueryFrame,blastReportDB-method}
\alias{getQueryFrame-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQueryFrame(x, max = FALSE)

  \S4method{getQueryFrame}{Hsp}(x, ...)

  \S4method{getQueryFrame}{HspList}(x, ...)

  \S4method{getQueryFrame}{Hit}(x, ...)

  \S4method{getQueryFrame}{HitList}(x, ...)

  \S4method{getQueryFrame}{Iteration}(x, ...)

  \S4method{getQueryFrame}{blastReportDB}(x, ...)
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
  Extract the frame of the query sequence.
}

