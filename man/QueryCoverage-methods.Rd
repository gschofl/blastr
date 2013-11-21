\docType{methods}
\name{getQueryCoverage}
\alias{getQueryCoverage}
\alias{getQueryCoverage,Hit-method}
\alias{getQueryCoverage,HitList-method}
\alias{getQueryCoverage,Hsp-method}
\alias{getQueryCoverage,HspList-method}
\alias{getQueryCoverage,Iteration-method}
\alias{getQueryCoverage,IterationList-method}
\alias{getQueryCoverage,blastReportDB-method}
\title{Extract the total Query Coverage.}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}
}
\description{
  Caluculated as (alignment length - gaps) / query length.
  If hsps overlap, they are merged before caluclating the
  coverage.
}

