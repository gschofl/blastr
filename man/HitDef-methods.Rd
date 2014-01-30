\docType{methods}
\name{getHitDef,Hit-method}
\alias{getHitDef}
\alias{getHitDef,Hit-method}
\alias{getHitDef,HitList-method}
\alias{getHitDef,Iteration-method}
\alias{getHitDef,blastReportDB-method}
\alias{getPrimaryHitDef}
\alias{getPrimaryHitDef,Hit-method}
\alias{getPrimaryHitDef,HitList-method}
\alias{getPrimaryHitDef,Iteration-method}
\alias{getPrimaryHitDef,IterationList-method}
\title{Access components of Blast Hits}
\usage{
\S4method{getHitDef}{Hit}(x)

\S4method{getHitDef}{HitList}(x)

\S4method{getPrimaryHitDef}{Hit}(x)

\S4method{getPrimaryHitDef}{HitList}(x)

\S4method{getHitDef}{Iteration}(x)

\S4method{getPrimaryHitDef}{Iteration}(x)

\S4method{getPrimaryHitDef}{IterationList}(x)

getHitDef(x, ...)

getPrimaryHitDef(x, ...)

\S4method{getHitDef}{blastReportDB}(x, id, ...)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
A character vector or list of character vectors.
}
\description{
Extract the (primary) hit definition.
}

