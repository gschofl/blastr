\docType{methods}
\name{getIdentity,Hit-method}
\alias{getIdentity}
\alias{getIdentity,Hit-method}
\alias{getIdentity,HitList-method}
\alias{getIdentity,Hsp-method}
\alias{getIdentity,HspList-method}
\alias{getIdentity,Iteration-method}
\alias{getIdentity,blastReportDB-method}
\alias{getPositive,Iteration-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getIdentity}{Hit}(x, max = FALSE)

\S4method{getIdentity}{HitList}(x, max = FALSE)

\S4method{getIdentity}{Hsp}(x, max = FALSE)

\S4method{getIdentity}{HspList}(x, max = FALSE)

\S4method{getIdentity}{Iteration}(x, max = FALSE)

\S4method{getPositive}{Iteration}(x, max = FALSE)

getIdentity(x, max = FALSE)

\S4method{getIdentity}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the numbers of identities.
}

