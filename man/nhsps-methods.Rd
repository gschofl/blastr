\docType{methods}
\name{nhsps,Hit-method}
\alias{nhsps}
\alias{nhsps,Hit-method}
\alias{nhsps,HitList-method}
\alias{nhsps,Iteration-method}
\alias{nhsps,IterationList-method}
\title{Access components of Blast Hits}
\usage{
\S4method{nhsps}{Hit}(x)

\S4method{nhsps}{HitList}(x)

\S4method{nhsps}{Iteration}(x)

\S4method{nhsps}{IterationList}(x)

nhsps(x)
}
\arguments{
  \item{x}{A \linkS4class{Hit} or \linkS4class{HitList}
  object.}
}
\value{
An numeric vector.
}
\description{
Extract the number of Hsps for each Hit.
}

