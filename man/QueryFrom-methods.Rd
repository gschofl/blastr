\name{getQueryFrom}
\alias{getQueryFrom}
\alias{getQueryFrom,Hit-method}
\alias{getQueryFrom,HitList-method}
\alias{getQueryFrom,Hsp-method}
\alias{getQueryFrom,HspList-method}
\alias{getQueryFrom,Iteration-method}
\alias{getQueryFrom,blastReportDB-method}
\alias{getQueryFrom-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQueryFrom(x, max = FALSE)

  \S4method{getQueryFrom}{Hsp}(x, ...)

  \S4method{getQueryFrom}{HspList}(x, ...)

  \S4method{getQueryFrom}{Hit}(x, ...)

  \S4method{getQueryFrom}{HitList}(x, ...)

  \S4method{getQueryFrom}{Iteration}(x, ...)

  \S4method{getQueryFrom}{blastReportDB}(x, ...)
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
  Extract the position of the start residues the query
  sequence.
}

