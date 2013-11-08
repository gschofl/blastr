\name{getEvalue}
\alias{getEvalue}
\alias{getEvalue,Hit-method}
\alias{getEvalue,HitList-method}
\alias{getEvalue,Hsp-method}
\alias{getEvalue,HspList-method}
\alias{getEvalue,Iteration-method}
\alias{getEvalue,blastReportDB-method}
\alias{getEvalue-methods}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
  getEvalue(x, max = FALSE)

  \S4method{getEvalue}{Hsp}(x, ...)

  \S4method{getEvalue}{HspList}(x, ...)

  \S4method{getEvalue}{Hit}(x, ...)

  \S4method{getEvalue}{HitList}(x, ...)

  \S4method{getEvalue}{Iteration}(x, ...)

  \S4method{getEvalue}{blastReportDB}(x, ...)
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

