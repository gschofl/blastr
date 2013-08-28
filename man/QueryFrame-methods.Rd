\docType{methods}
\name{getQueryFrame}
\alias{getQueryFrame}
\alias{getQueryFrame,Hit-method}
\alias{getQueryFrame,HitList-method}
\alias{getQueryFrame,Hsp-method}
\alias{getQueryFrame,HspList-method}
\alias{getQueryFrame,Iteration-method}
\alias{getQueryFrame,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQueryFrame(x, max = FALSE)
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

