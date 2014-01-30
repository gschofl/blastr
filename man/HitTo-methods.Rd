\docType{methods}
\name{getHitTo,Hit-method}
\alias{getHitTo}
\alias{getHitTo,Hit-method}
\alias{getHitTo,HitList-method}
\alias{getHitTo,Hsp-method}
\alias{getHitTo,HspList-method}
\alias{getHitTo,Iteration-method}
\alias{getHitTo,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getHitTo}{Hit}(x, max = FALSE)

\S4method{getHitTo}{HitList}(x, max = FALSE)

\S4method{getHitTo}{Hsp}(x, max = FALSE)

\S4method{getHitTo}{HspList}(x, max = FALSE)

\S4method{getHitTo}{Iteration}(x, max = FALSE)

getHitTo(x, max = FALSE)

\S4method{getHitTo}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the position of end residues for the hit sequence.
}

