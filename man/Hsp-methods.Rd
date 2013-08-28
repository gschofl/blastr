\docType{methods}
\name{getHsp}
\alias{getHsp}
\alias{getHsp,Hit-method}
\alias{getHsp,HitList-method}
\alias{getHsp,Iteration-method}
\title{Access components of Blast Hits}
\usage{
  getHsp(x, ...)
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

