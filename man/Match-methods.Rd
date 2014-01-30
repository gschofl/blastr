\docType{methods}
\name{getMatch,Hit-method}
\alias{getMatch}
\alias{getMatch,Hit-method}
\alias{getMatch,HitList-method}
\alias{getMatch,Hsp-method}
\alias{getMatch,HspList-method}
\alias{getMatch,Iteration-method}
\alias{getMatch,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getMatch}{Hit}(x, max = FALSE)

\S4method{getMatch}{HitList}(x, max = FALSE)

\S4method{getMatch}{Hsp}(x, max = FALSE)

\S4method{getMatch}{HspList}(x, max = FALSE)

\S4method{getMatch}{Iteration}(x, max = FALSE)

getMatch(x, max = FALSE)

\S4method{getMatch}{blastReportDB}(x, id, max = FALSE, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
An \linkS4class{BString} or \linkS4class{BStringSet}
object.
}
\description{
Extract the match sequence (midline).
}

