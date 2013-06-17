\name{getGeneID}
\alias{getGeneID}
\alias{getGeneID,Hit-method}
\alias{getGeneID,HitList-method}
\alias{getGeneID,Iteration-method}
\alias{getGeneID,IterationList-method}
\alias{getGeneID,blastReportDB-method}
\alias{getGeneID-methods}
\title{Access components of Blast Hits}
\usage{
  getGeneID(x, ...)

  \S4method{getGeneID}{Hit}(x, ...)

  \S4method{getGeneID}{HitList}(x, ...)

  \S4method{getGeneID}{Iteration}(x, ...)

  \S4method{getGeneID}{IterationList}(x, ...)

  \S4method{getGeneID}{blastReportDB}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  A character vector or list of character vectors.
}
\description{
  Extract the GIs (gene identifiers) of hits.
}

