\docType{methods}
\name{getHitRange}
\alias{getHitRange}
\alias{getHitRange,Hit-method}
\alias{getHitRange,HitList-method}
\alias{getHitRange,Hsp-method}
\alias{getHitRange,HspList-method}
\alias{getHitRange,Iteration-method}
\alias{getHitRange,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
  An \linkS4class{IRanges} object.
}
\description{
  Extract the alignment ranges for the hit sequence.
}

