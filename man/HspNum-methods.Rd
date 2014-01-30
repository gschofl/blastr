\docType{methods}
\name{getHspNum,Hit-method}
\alias{getHspNum}
\alias{getHspNum,Hit-method}
\alias{getHspNum,HitList-method}
\alias{getHspNum,Hsp-method}
\alias{getHspNum,HspList-method}
\alias{getHspNum,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getHspNum}{Hit}(x, max = FALSE)

\S4method{getHspNum}{HitList}(x, max = FALSE)

\S4method{getHspNum}{Hsp}(x, max = FALSE)

\S4method{getHspNum}{HspList}(x, max = FALSE)

getHspNum(x, ...)

getHspNum(x, id)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return index of the Hsp with
  the highest bit score.}
}
\value{
An numeric vector or a list of numeric vectors
}
\description{
Extract the Hsp indices.
}

