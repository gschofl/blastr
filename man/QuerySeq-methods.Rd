\docType{methods}
\name{getQuerySeq}
\alias{getQuerySeq}
\alias{getQuerySeq,Hit-method}
\alias{getQuerySeq,HitList-method}
\alias{getQuerySeq,Hsp-method}
\alias{getQuerySeq,HspList-method}
\alias{getQuerySeq,Iteration-method}
\alias{getQuerySeq,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getQuerySeq(x, max = FALSE)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
  An \linkS4class{BString} or \linkS4class{BStringSet}
  object.
}
\description{
  Extract the query sequence.
}

