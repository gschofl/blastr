\name{getIdentity}
\alias{getIdentity}
\alias{getIdentity,Hit-method}
\alias{getIdentity,HitList-method}
\alias{getIdentity,Hsp-method}
\alias{getIdentity,HspList-method}
\alias{getIdentity,Iteration-method}
\alias{getIdentity,blastReportDB-method}
\alias{getIdentity-methods}
\alias{getPositive,Iteration-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getIdentity(x, max = FALSE)

  \S4method{getIdentity}{Hsp}(x, ...)

  \S4method{getIdentity}{HspList}(x, ...)

  \S4method{getIdentity}{Hit}(x, ...)

  \S4method{getIdentity}{HitList}(x, ...)

  \S4method{getIdentity}{Iteration}(x, ...)

  \S4method{getPositive}{Iteration}(x, ...)

  \S4method{getIdentity}{blastReportDB}(x, ...)
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

