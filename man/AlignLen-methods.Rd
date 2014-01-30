\docType{methods}
\name{getAlignLen,Hit-method}
\alias{getAlignLen}
\alias{getAlignLen,Hit-method}
\alias{getAlignLen,HitList-method}
\alias{getAlignLen,Hsp-method}
\alias{getAlignLen,HspList-method}
\alias{getAlignLen,Iteration-method}
\alias{getAlignLen,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getAlignLen}{Hit}(x, max = FALSE)

\S4method{getAlignLen}{HitList}(x, max = FALSE)

\S4method{getAlignLen}{Hsp}(x, max = FALSE)

\S4method{getAlignLen}{HspList}(x, max = FALSE)

\S4method{getAlignLen}{Iteration}(x, max = FALSE)

getIdentity(x, max = FALSE)

\S4method{getAlignLen}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the alignment lengths.
}

