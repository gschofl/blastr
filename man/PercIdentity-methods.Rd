\name{getPercIdentity}
\alias{getMaxPercIdentity}
\alias{getMaxPercIdentity,Hit-method}
\alias{getMaxPercIdentity,HitList-method}
\alias{getMaxPercIdentity,Hsp-method}
\alias{getMaxPercIdentity,HspList-method}
\alias{getMaxPercIdentity,Iteration-method}
\alias{getMaxPercIdentity,IterationList-method}
\alias{getMaxPercIdentity,blastReportDB-method}
\alias{getMaxPercIdentity-methods}
\alias{getPercIdentity}
\alias{getPercIdentity,Hit-method}
\alias{getPercIdentity,HitList-method}
\alias{getPercIdentity,Hsp-method}
\alias{getPercIdentity,HspList-method}
\alias{getPercIdentity,Iteration-method}
\alias{getPercIdentity,blastReportDB-method}
\alias{getPercIdentity-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getPercIdentity(x, max = FALSE)

  getMaxPercIdentity(x)

  \S4method{getPercIdentity}{Hsp}(x, ...)

  \S4method{getPercIdentity}{HspList}(x, ...)

  \S4method{getMaxPercIdentity}{Hsp}(x, ...)

  \S4method{getMaxPercIdentity}{HspList}(x, ...)

  \S4method{getPercIdentity}{Hit}(x, ...)

  \S4method{getPercIdentity}{HitList}(x, ...)

  \S4method{getMaxPercIdentity}{Hit}(x, ...)

  \S4method{getMaxPercIdentity}{HitList}(x, ...)

  \S4method{getPercIdentity}{Iteration}(x, ...)

  \S4method{getMaxPercIdentity}{Iteration}(x, ...)

  \S4method{getMaxPercIdentity}{IterationList}(x, ...)

  \S4method{getPercIdentity}{blastReportDB}(x, ...)

  \S4method{getMaxPercIdentity}{blastReportDB}(x, ...)
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
  Extract Percent Identity, Percent Positives or Percent
  Gaps.
}

