\docType{methods}
\name{getHitCoverage}
\alias{getHitCoverage}
\alias{getHitCoverage,Hit-method}
\alias{getHitCoverage,HitList-method}
\alias{getHitCoverage,Hsp-method}
\alias{getHitCoverage,HspList-method}
\alias{getHitCoverage,Iteration-method}
\alias{getHitCoverage,IterationList-method}
\alias{getHitCoverage,blastReportDB-method}
\title{Extract the total Hit Coverage.}
\usage{
  getHitCoverage(x)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}
}
\description{
  Caluculated as (alignment length - gaps) / hit length. If
  hsps overlap, they are merged before caluclating the
  coverage.
}

