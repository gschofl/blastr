\docType{methods}
\name{getQuerySeq,Hit-method}
\alias{getQuerySeq}
\alias{getQuerySeq,Hit-method}
\alias{getQuerySeq,HitList-method}
\alias{getQuerySeq,Hsp-method}
\alias{getQuerySeq,HspList-method}
\alias{getQuerySeq,Iteration-method}
\alias{getQuerySeq,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getQuerySeq}{Hit}(x, max = FALSE)

\S4method{getQuerySeq}{HitList}(x, max = FALSE)

\S4method{getQuerySeq}{Hsp}(x, max = FALSE)

\S4method{getQuerySeq}{HspList}(x, max = FALSE)

\S4method{getQuerySeq}{Iteration}(x, max = FALSE)

getQuerySeq(x, max = FALSE)

\S4method{getQuerySeq}{blastReportDB}(x, id, max = FALSE, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}

  \item{max}{If \code{TRUE} return data only for the hsp
  with the highest bit score.}
}
\value{
An \linkS4class{BString} or \linkS4class{BStringSet}
object.
}
\description{
Extract the query sequence.
}

