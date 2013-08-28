\docType{methods}
\name{getPercIdentity}
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
\alias{getPercIdentity,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getPercIdentity(x, max = FALSE)

  getMaxPercIdentity(x)
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

