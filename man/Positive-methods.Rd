\docType{methods}
\name{getPositive,Hit-method}
\alias{getPositive}
\alias{getPositive,Hit-method}
\alias{getPositive,HitList-method}
\alias{getPositive,Hsp-method}
\alias{getPositive,HspList-method}
\alias{getPositive,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getPositive}{Hit}(x, max = FALSE)

\S4method{getPositive}{HitList}(x, max = FALSE)

\S4method{getPositive}{Hsp}(x, max = FALSE)

\S4method{getPositive}{HspList}(x, max = FALSE)

getPositive(x, max = FALSE)

\S4method{getPositive}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the numbers of positives.
}

