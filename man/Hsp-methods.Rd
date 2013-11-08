\name{getHsp}
\alias{getHsp}
\alias{getHsp,Hit-method}
\alias{getHsp,HitList-method}
\alias{getHsp,Iteration-method}
\alias{getHsp-methods}
\title{Access components of Blast Hits}
\usage{
  getHsp(x, ...)

  \S4method{getHsp}{Hit}(x, ...)

  \S4method{getHsp}{HitList}(x, ...)

  \S4method{getHsp}{Iteration}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
  An \linkS4class{Hsp} or \linkS4class{HspList}.
}
\description{
  Extract Hsps from a Hit.
}

