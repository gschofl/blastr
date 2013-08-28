\docType{methods}
\name{getScore}
\alias{getScore}
\alias{getScore,Hit-method}
\alias{getScore,HitList-method}
\alias{getScore,Hsp-method}
\alias{getScore,HspList-method}
\alias{getScore,Iteration-method}
\alias{getScore,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getScore(x, max = FALSE)
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

