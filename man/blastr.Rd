\docType{package}
\name{blastr}
\alias{blastr}
\alias{blastr-package}
\title{A wrapper around NCBI Blast.}
\description{
Provide a set of classes and functions to perform and parse
BLAST queries locally or vie NCBI Blast web service
queries.
}
\details{
\tabular{ll}{ Package: \tab blastr\cr Type: \tab Package\cr
Version: \tab 0.2.1\cr Date: \tab 2013-11-22\cr License:
\tab MIT\cr }
}
\section{Package options}{
  \emph{blastr} uses the following \code{\link{options}} to
  configure behaviour:

  \itemize{ \item \code{blastr.blastdb.path}: Path to a
  local installation of the NCBI BLAST databases. By
  default this option is not set. I recommend setting it to
  a path of your choice in your .Rprofile file. Run
  \code{\link{update_blastdb}} to install these databases
  locally. }
}
\author{
Gerhard Schoefl \email{gschofl@yahoo.de}

Maintainer: Gerhard Schoefl \email{gschofl@yahoo.de}
}
\keyword{package}

