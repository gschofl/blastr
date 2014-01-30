\docType{methods}
\name{getQueryTo,Hit-method}
\alias{getQueryTo}
\alias{getQueryTo,Hit-method}
\alias{getQueryTo,HitList-method}
\alias{getQueryTo,Hsp-method}
\alias{getQueryTo,HspList-method}
\alias{getQueryTo,Iteration-method}
\alias{getQueryTo,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getQueryTo}{Hit}(x, max = FALSE)

\S4method{getQueryTo}{HitList}(x, max = FALSE)

\S4method{getQueryTo}{Hsp}(x, max = FALSE)

\S4method{getQueryTo}{HspList}(x, max = FALSE)

\S4method{getQueryTo}{Iteration}(x, max = FALSE)

getQueryTo(x, max = FALSE)

\S4method{getQueryTo}{blastReportDB}(x, id, max = FALSE, ...)
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

