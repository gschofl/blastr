\name{getAccession}
\alias{getAccession}
\alias{getAccession,Hit-method}
\alias{getAccession,HitList-method}
\alias{getAccession,Iteration-method}
\alias{getAccession,IterationList-method}
\alias{getAccession,blastReportDB-method}
\alias{getAccession-methods}
\title{Access components of Blast Hits}
\usage{
  getAccession(x, ...)

  \S4method{getAccession}{Hit}(x, ...)

  \S4method{getAccession}{HitList}(x, ...)

  \S4method{getAccession}{Iteration}(x, ...)

  \S4method{getAccession}{IterationList}(x, ...)

  \S4method{getAccession}{blastReportDB}(x, ...)
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
  Extract the accession numbers of hits.
}

