\docType{methods}
\name{getHitNum}
\alias{HitNum,Iteration-method}
\alias{getHitNum}
\alias{getHitNum,Hit-method}
\alias{getHitNum,HitList-method}
\alias{getHitNum,IterationList-method}
\alias{getHitNum,blastReportDB-method}
\title{Access components of Blast Hits}
\usage{
  getHitNum(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  An \linkS4class{Hsp} or \linkS4class{HspList}.

  A numeric vector or list of numeric vectors.
}
\description{
  Extract the Hit indices.
}

