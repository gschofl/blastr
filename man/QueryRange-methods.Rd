\docType{methods}
\name{getQueryRange}
\alias{getQueryRange}
\alias{getQueryRange,Hit-method}
\alias{getQueryRange,Hsp-method}
\alias{getQueryRange,HspList-method}
\alias{getQueryRange,Iteration-method}
\alias{getQueryRange,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQueryRange(x, max = FALSE)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
  An \linkS4class{IRanges} object.
}
\description{
  Extract the alignment ranges for the query sequence.
}

