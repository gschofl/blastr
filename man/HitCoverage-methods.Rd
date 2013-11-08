\name{getHitCoverage}
\alias{getHitCoverage}
\alias{getHitCoverage,Hit-method}
\alias{getHitCoverage,HitList-method}
\alias{getHitCoverage,Hsp-method}
\alias{getHitCoverage,HspList-method}
\alias{getHitCoverage,Iteration-method}
\alias{getHitCoverage,IterationList-method}
\alias{getHitCoverage,blastReportDB-method}
\alias{getHitCoverage-methods}
\title{Extract the total Hit Coverage.}
\usage{
  getHitCoverage(x)

  \S4method{getHitCoverage}{Hsp}(x, ...)

  \S4method{getHitCoverage}{HspList}(x, ...)

  \S4method{getHitCoverage}{Hit}(x, ...)

  \S4method{getHitCoverage}{HitList}(x, ...)

  \S4method{getHitCoverage}{Iteration}(x, ...)

  \S4method{getHitCoverage}{IterationList}(x, ...)

  \S4method{getHitCoverage}{blastReportDB}(x, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}
}
\description{
  Caluculated as (alignment length - gaps) / hit length. If
  hsps overlap, they are merged before caluclating the
  coverage.
}

