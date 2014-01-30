\docType{methods}
\name{getHitNum,Hit-method}
\alias{HitNum,Iteration-method}
\alias{getHitNum}
\alias{getHitNum,Hit-method}
\alias{getHitNum,HitList-method}
\alias{getHitNum,Iteration-method}
\alias{getHitNum,IterationList-method}
\alias{getHitNum,blastReportDB-method}
\title{Access components of Blast Hits}
\usage{
\S4method{getHitNum}{Hit}(x)

\S4method{getHitNum}{HitList}(x)

\S4method{getHitNum}{Iteration}(x)

\S4method{getHitNum}{IterationList}(x)

getHitNum(x, ...)

\S4method{getHitNum}{blastReportDB}(x, id, ...)
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

