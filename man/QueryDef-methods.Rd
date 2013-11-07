\name{getQueryDef}
\alias{getQueryDef}
\alias{getQueryDef,Iteration-method}
\alias{getQueryDef,IterationList-method}
\alias{getQueryDef,blastReport-method}
\alias{getQueryDef,blastReportDB-method}
\alias{getQueryDef,blastReportDb-method}
\alias{getQueryDef-methods}
\title{Access components of Blast Iterations (Queries)}
\usage{
  getQueryDef(x, ...)

  \S4method{getQueryDef}{Iteration}(x, ...)

  \S4method{getQueryDef}{IterationList}(x, ...)

  \S4method{getQueryDef}{blastReport}(x, ...)

  \S4method{getQueryDef}{blastReportDB}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{blastReport} or
  \linkS4class{IterationList} object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  A character vector or list of character vectors.
}
\description{
  Extract the query definition.
}

