\docType{methods}
\name{getBitscore}
\alias{getBitscore}
\alias{getBitscore,Hit-method}
\alias{getBitscore,HitList-method}
\alias{getBitscore,Hsp-method}
\alias{getBitscore,HspList-method}
\alias{getBitscore,Iteration-method}
\alias{getBitscore,blastReportDB-method}
\alias{getMaxBitscore}
\alias{getMaxBitscore,Hit-method}
\alias{getMaxBitscore,HitList-method}
\alias{getMaxBitscore,Iteration-method}
\alias{getMaxBitscore,IterationList-method}
\alias{getMaxBitscore,blastReportDB-method}
\alias{getTotalBitscore}
\alias{getTotalBitscore,Hit-method}
\alias{getTotalBitscore,HitList-method}
\alias{getTotalBitscore,Iteration-method}
\alias{getTotalBitscore,IterationList-method}
\alias{getTotalBitscore,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getBitscore(x, max = FALSE, sum = FALSE)

  getMaxBitscore(x)

  getTotalBitscore(x)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return only the highest score.}

  \item{sum}{If \code{TRUE} sum scores from all hsps of a
  Blast hit.}
}
\value{
  An numeric vector or a list of numeric vectors
}
\description{
  Extract the bit score.
}

