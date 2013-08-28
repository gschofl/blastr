\docType{methods}
\name{show}
\alias{show}
\alias{show,Hit-method}
\alias{show,HitList-method}
\alias{show,Hsp-method}
\alias{show,HspList-method}
\alias{show,Iteration-method}
\alias{show,IterationList-method}
\alias{show,blastReport-method}
\alias{show,blastReportDB-method}
\alias{show,blastTable-method}
\title{Method extensions to show for blast objects.}
\description{
  The \code{show} methods for various blast objects can be
  modified by a number of global options. Specifically, the
  number of \code{\linkS4class{Hit}s} shown when displaying
  \code{\linkS4class{Iteration}s} or
  \code{\linkS4class{HitList}s} are set by \code{showHits}.
  Whether alignments are displayed or not is controlled by
  \code{showAlignment}. Setting these options to
  \code{NULL}, restores the defaults.
}
\details{
  See the general documentation of
  \code{\link[methods]{show}} method for the expected
  behavior.
}
\examples{
# options("showHits" = 20)
# show(hitlist)
# options("showHits" = NULL)
# show(hitlist)
}
\seealso{
  \code{\link[methods]{show}}
}

