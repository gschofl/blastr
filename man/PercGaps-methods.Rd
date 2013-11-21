\docType{methods}
\name{getPercGaps}
\alias{getPercGaps}
\alias{getPercGaps,Hit-method}
\alias{getPercGaps,HitList-method}
\alias{getPercGaps,Hsp-method}
\alias{getPercGaps,HspList-method}
\alias{getPercGaps,Iteration-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
  An numeric vector or a list of numeric vectors
}
\description{
  Extract Percent Gaps.
}

