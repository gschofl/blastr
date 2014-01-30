\name{blastReport}
\alias{blastReport}
\title{Parse NCBI BLAST XML files into \linkS4class{blastReport} objects.}
\usage{
blastReport(blast, asText = FALSE)
}
\arguments{
  \item{blastfile}{Blast output in XML format (file path or
  character vector)}

  \item{asText}{If \code{TRUE} the XML blast output is
  passed in as a string instead of a file path.}
}
\value{
A \code{\linkS4class{blastReport}} object.
}
\description{
Parse NCBI BLAST XML files into \linkS4class{blastReport}
objects.
}

