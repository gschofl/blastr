\docType{methods}
\name{getHitLen,Hit-method}
\alias{getHitLen}
\alias{getHitLen,Hit-method}
\alias{getHitLen,HitList-method}
\alias{getHitLen,Iteration-method}
\alias{getHitLen,IterationList-method}
\alias{getHitLen,blastReportDB-method}
\title{Access components of Blast Hits}
\usage{
\S4method{getHitLen}{Hit}(x)

\S4method{getHitLen}{HitList}(x)

\S4method{getHitLen}{Iteration}(x)

\S4method{getHitLen}{IterationList}(x)

getHitLen(x, ...)

\S4method{getHitLen}{blastReportDB}(x, id, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
A numeric vector or list of numeric vectors.
}
\description{
Extract the hit length.
}

