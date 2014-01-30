\docType{methods}
\name{getGaps,Hit-method}
\alias{getGaps}
\alias{getGaps,Hit-method}
\alias{getGaps,HitList-method}
\alias{getGaps,Hsp-method}
\alias{getGaps,HspList-method}
\alias{getGaps,Iteration-method}
\alias{getGaps,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getGaps}{Hit}(x, max = FALSE)

\S4method{getGaps}{HitList}(x, max = FALSE)

\S4method{getGaps}{Hsp}(x, max = FALSE)

\S4method{getGaps}{HspList}(x, max = FALSE)

\S4method{getGaps}{Iteration}(x, max = FALSE)

getGaps(x, max = FALSE)

\S4method{getGaps}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the numbers of gaps.
}

