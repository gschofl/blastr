\docType{methods}
\name{getScore,Hit-method}
\alias{getScore}
\alias{getScore,Hit-method}
\alias{getScore,HitList-method}
\alias{getScore,Hsp-method}
\alias{getScore,HspList-method}
\alias{getScore,Iteration-method}
\alias{getScore,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getScore}{Hit}(x, max = FALSE)

\S4method{getScore}{HitList}(x, max = FALSE)

\S4method{getScore}{Hsp}(x, max = FALSE)

\S4method{getScore}{HspList}(x, max = FALSE)

\S4method{getScore}{Iteration}(x, max = FALSE)

getScore(x, max = FALSE)

\S4method{getScore}{blastReportDB}(x, id, max = FALSE, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return only the highest score.}
}
\value{
An numeric vector or a list of numeric vectors
}
\description{
Extract the raw score.
}

