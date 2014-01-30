\docType{methods}
\name{getPercPositive,Hit-method}
\alias{getPercPositive}
\alias{getPercPositive,Hit-method}
\alias{getPercPositive,HitList-method}
\alias{getPercPositive,Hsp-method}
\alias{getPercPositive,HspList-method}
\alias{getPercPositive,Iteration-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getPercPositive}{Hit}(x, max = FALSE)

\S4method{getPercPositive}{HitList}(x, max = FALSE)

\S4method{getPercPositive}{Hsp}(x, max = FALSE)

\S4method{getPercPositive}{HspList}(x, max = FALSE)

\S4method{getPercPositive}{Iteration}(x, max = FALSE)

getPercPositive(x, max = FALSE)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
An numeric vector or a list of numeric vectors
}
\description{
Extract Percent Positives.
}

