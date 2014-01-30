\docType{methods}
\name{getQueryCoverage,Hit-method}
\alias{getQueryCoverage}
\alias{getQueryCoverage,Hit-method}
\alias{getQueryCoverage,HitList-method}
\alias{getQueryCoverage,Hsp-method}
\alias{getQueryCoverage,HspList-method}
\alias{getQueryCoverage,Iteration-method}
\alias{getQueryCoverage,IterationList-method}
\alias{getQueryCoverage,blastReportDB-method}
\title{Extract the total Query Coverage.}
\usage{
\S4method{getQueryCoverage}{Hit}(x)

\S4method{getQueryCoverage}{HitList}(x)

\S4method{getQueryCoverage}{Hsp}(x)

\S4method{getQueryCoverage}{HspList}(x)

\S4method{getQueryCoverage}{Iteration}(x)

\S4method{getQueryCoverage}{IterationList}(x)

getQueryCoverage(x)

\S4method{getQueryCoverage}{blastReportDB}(x, id, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hsp} or \linkS4class{HspList}.}
}
\description{
Caluculated as (alignment length - gaps) / query length. If
hsps overlap, they are merged before caluclating the
coverage.
}

