\name{blastReport}
\alias{blastReport}
\title{Parse XML NCBI BLAST files}
\usage{
  blastReport(blast_file, asText = FALSE)
}
\arguments{
  \item{blast_file}{Blast output in XML format (File or
  character vector)}

  \item{asText}{If \code{TRUE} the XML blast output is
  passed in as a string instead of a file path.}
}
\value{
  A \code{\linkS4class{blastReport}} object.
}
\description{
  Parse XML NCBI BLAST files
}

