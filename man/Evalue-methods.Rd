\docType{methods}
\name{getEvalue}
\alias{getEvalue}
\alias{getEvalue,Hit-method}
\alias{getEvalue,HitList-method}
\alias{getEvalue,Hsp-method}
\alias{getEvalue,HspList-method}
\alias{getEvalue,Iteration-method}
\alias{getEvalue,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getEvalue(x, max = FALSE)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return only the highest
  evalue.}
}
\value{
  An numeric vector or a list of numeric vectors
}
\description{
  Extract the expect value.
}

