\name{getQueryTo}
\alias{getQueryTo}
\alias{getQueryTo,Hit-method}
\alias{getQueryTo,HitList-method}
\alias{getQueryTo,Hsp-method}
\alias{getQueryTo,HspList-method}
\alias{getQueryTo,Iteration-method}
\alias{getQueryTo,blastReportDB-method}
\alias{getQueryTo-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQueryTo(x, max = FALSE)

  \S4method{getQueryTo}{Hsp}(x, ...)

  \S4method{getQueryTo}{HspList}(x, ...)

  \S4method{getQueryTo}{Hit}(x, ...)

  \S4method{getQueryTo}{HitList}(x, ...)

  \S4method{getQueryTo}{Iteration}(x, ...)

  \S4method{getQueryTo}{blastReportDB}(x, ...)
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
  Extract the position of end residues for the query
  sequence.
}

