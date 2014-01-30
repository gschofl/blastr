\docType{methods}
\name{getPercIdentity,Hit-method}
\alias{getMaxPercIdentity}
\alias{getMaxPercIdentity,Hit-method}
\alias{getMaxPercIdentity,HitList-method}
\alias{getMaxPercIdentity,Hsp-method}
\alias{getMaxPercIdentity,HspList-method}
\alias{getMaxPercIdentity,Iteration-method}
\alias{getMaxPercIdentity,IterationList-method}
\alias{getMaxPercIdentity,blastReportDB-method}
\alias{getPercIdentity}
\alias{getPercIdentity,Hit-method}
\alias{getPercIdentity,HitList-method}
\alias{getPercIdentity,Hsp-method}
\alias{getPercIdentity,HspList-method}
\alias{getPercIdentity,Iteration-method}
\alias{getPercIdentity,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getPercIdentity}{Hit}(x, max = FALSE)

\S4method{getPercIdentity}{HitList}(x, max = FALSE)

\S4method{getMaxPercIdentity}{Hit}(x)

\S4method{getMaxPercIdentity}{HitList}(x)

\S4method{getPercIdentity}{Hsp}(x, max = FALSE)

\S4method{getPercIdentity}{HspList}(x, max = FALSE)

\S4method{getMaxPercIdentity}{Hsp}(x)

\S4method{getMaxPercIdentity}{HspList}(x)

\S4method{getPercIdentity}{Iteration}(x, max = FALSE)

\S4method{getMaxPercIdentity}{Iteration}(x)

\S4method{getMaxPercIdentity}{IterationList}(x)

getPercIdentity(x, max = FALSE)

getMaxPercIdentity(x)

\S4method{getPercIdentity}{blastReportDB}(x, id, max = FALSE, ...)

\S4method{getMaxPercIdentity}{blastReportDB}(x, id, ...)
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

