\docType{methods}
\name{getQueryTo}
\alias{getQueryTo}
\alias{getQueryTo,Hit-method}
\alias{getQueryTo,HitList-method}
\alias{getQueryTo,Hsp-method}
\alias{getQueryTo,HspList-method}
\alias{getQueryTo,Iteration-method}
\alias{getQueryTo,blastReportDB-method}
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
  Extract the position of end residues for the query
  sequence.
}

