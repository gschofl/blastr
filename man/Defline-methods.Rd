\name{getDefline}
\alias{getDefline}
\alias{getDefline,Hit-method}
\alias{getDefline,HitList-method}
\alias{getDefline,Iteration-method}
\alias{getDefline-methods}
\alias{getPrimaryDefline}
\alias{getPrimaryDefline,Hit-method}
\alias{getPrimaryDefline,HitList-method}
\alias{getPrimaryDefline,Iteration-method}
\alias{getPrimaryDefline,IterationList-method}
\alias{getPrimaryDefline-methods}
\title{Access components of Blast Hits}
\usage{
  getDefline(x, ...)

  getPrimaryDefline(x, ...)

  \S4method{getDefline}{Hit}(x, ...)

  \S4method{getDefline}{HitList}(x, ...)

  \S4method{getPrimaryDefline}{Hit}(x, ...)

  \S4method{getPrimaryDefline}{HitList}(x, ...)

  \S4method{getDefline}{Iteration}(x, ...)

  \S4method{getPrimaryDefline}{Iteration}(x, ...)

  \S4method{getPrimaryDefline}{IterationList}(x, ...)
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
  Extract the (primary) definition lines.
}

