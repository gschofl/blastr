\docType{methods}
\name{getQueryFrame,Hit-method}
\alias{getQueryFrame}
\alias{getQueryFrame,Hit-method}
\alias{getQueryFrame,HitList-method}
\alias{getQueryFrame,Hsp-method}
\alias{getQueryFrame,HspList-method}
\alias{getQueryFrame,Iteration-method}
\alias{getQueryFrame,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getQueryFrame}{Hit}(x, max = FALSE)

\S4method{getQueryFrame}{HitList}(x, max = FALSE)

\S4method{getQueryFrame}{Hsp}(x, max = FALSE)

\S4method{getQueryFrame}{HspList}(x, max = FALSE)

\S4method{getQueryFrame}{Iteration}(x, max = FALSE)

getQueryFrame(x, max = FALSE)

\S4method{getQueryFrame}{blastReportDB}(x, id, max = FALSE, ...)
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

