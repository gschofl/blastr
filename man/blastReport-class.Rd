\name{blastReport-class}
\alias{blastReport-class}
\title{blastReport-class}
\usage{
  new_blastReport(...)
}
\description{
  \code{"blastReport"} is the top-level container for data
  parsed from NCBI Blast XML output. It contains the
  component data classes: \code{\linkS4class{BlastHeader}},
  \code{\linkS4class{BlastParameters}}, and
  \code{\linkS4class{IterationList}}.
}
\details{
  The \code{\linkS4class{Iteration}} elements store results
  from individual BLAST queries. Each \code{Iteration}
  holds a \code{\linkS4class{HitList}} with possibly
  multiple \code{\linkS4class{Hit}s}, which, in turn, can
  contain multiple high-scoring pairs
  (\code{\linkS4class{Hsp}s}).

  Iterations, Hits, and Hsps can be extracted using the
  accessors \code{\link{getIteration}},
  \code{\link{getHit}}, and \code{\link{getHsp}}, or by
  directly subsetting a \code{blastReport} object.

  E.g. \code{report[[1]][[1]]} will return the first hit in
  the first query.
}
\section{Slots}{
  \describe{ \item{\code{header}:}{Header information;
  \code{"BlastHeader"}.} \item{\code{params}:}{Blast
  parameters and statistics \code{"BlastParameters"}.}
  \item{\code{iterations}:}{Iterations;
  \code{"IterationListOrChar"}.} }
}
\seealso{
  The constructor \code{\link{blastReport}}; the BLAST
  classes \code{\linkS4class{blastReportDB}} and
  \code{\linkS4class{blastTable}}
}

