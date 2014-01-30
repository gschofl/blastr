\docType{methods}
\name{getHsp,Hit-method}
\alias{getHsp}
\alias{getHsp,Hit-method}
\alias{getHsp,HitList-method}
\alias{getHsp,Iteration-method}
\title{Access components of Blast Hits}
\usage{
\S4method{getHsp}{Hit}(x, n = NULL, drop = TRUE)

\S4method{getHsp}{HitList}(x, n = NULL, drop = TRUE)

\S4method{getHsp}{Iteration}(x, n = NULL, drop = TRUE)

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

