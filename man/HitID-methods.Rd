\name{getHitID}
\alias{getHitID}
\alias{getHitID,Hit-method}
\alias{getHitID,HitList-method}
\alias{getHitID,Iteration-method}
\alias{getHitID,blastReportDB-method}
\alias{getHitID-methods}
\title{Access components of Blast Hits}
\usage{
  getHitID(x, db = 'any', ...)

  \S4method{getHitID}{Hit}(x, ...)

  \S4method{getHitID}{HitList}(x, ...)

  \S4method{getHitID}{Iteration}(x, ...)

  \S4method{getHitID}{blastReportDB}(x, ...)
}
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

