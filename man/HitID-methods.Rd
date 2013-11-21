\docType{methods}
\name{getHitID}
\alias{getHitID}
\alias{getHitID,Hit-method}
\alias{getHitID,HitList-method}
\alias{getHitID,Iteration-method}
\alias{getHitID,blastReportDB-method}
\title{Access components of Blast Hits}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{db}{Database tag (e.g. 'gi', 'ref', 'gb', 'dbj',
  'any')}

  \item{...}{Further arguments passed to methods.}
}
\value{
  An character matrix or character vector.
}
\description{
  Extract hit IDs.
}

