\docType{methods}
\name{getHitSeq,Hit-method}
\alias{getHitSeq}
\alias{getHitSeq,Hit-method}
\alias{getHitSeq,HitList-method}
\alias{getHitSeq,Hsp-method}
\alias{getHitSeq,HspList-method}
\alias{getHitSeq,Iteration-method}
\alias{getHitSeq,blastReportDB-method}
\title{Access components of High-Scoring Pairs (Hsps)}
\usage{
\S4method{getHitSeq}{Hit}(x, max = FALSE)

\S4method{getHitSeq}{HitList}(x, max = FALSE)

\S4method{getHitSeq}{Hsp}(x, max = FALSE)

\S4method{getHitSeq}{HspList}(x, max = FALSE)

\S4method{getHitSeq}{Iteration}(x, max = FALSE)

getHitSeq(x, max = FALSE)

\S4method{getHitSeq}{blastReportDB}(x, id, max = FALSE, ...)
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
Extract the hit sequence.
}

