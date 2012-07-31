\name{blast}
\alias{blast}
\title{Wrapper for the new NCBI BLAST+ tools}
\usage{
  blast(program = c("blastn", "blastp", "blastx", "tblastn", "tblastx"),
    query, db, outfmt, max_hits,
    strand = c("both", "plus", "minus"), ...,
    intern = FALSE, input = NULL, show_cmd = FALSE)
}
\arguments{
  \item{program}{One of blastn, blastp, blastx, tblastn, or
  tblastx}

  \item{query}{Query provided as an input file name, an
  \code{\link[Biostrings]{XString}} or
  \code{\link[Biostrings]{XStringSet}} object, or as a
  character vector.}

  \item{db}{Blast database name (defaults to 'nt' for
  blastn and 'nr' for the rest).}

  \item{outfmt}{Output format (defaults to XML output).}

  \item{max_hits}{Maximum number of hits to return.
  Defaults to 20. Sets the '-max_target_seqs' parameter
  internally.}

  \item{strand}{Query strand(s) to seach against database.}

  \item{...}{Arguments passed on to the blast commmand line
  tools.}

  \item{intern}{Set \code{TRUE} if no '-out' argument is
  specified. Captures the blast output in an R character
  vector.}

  \item{input}{Used to pass a character vector to the
  standard input.}

  \item{show_cmd}{If \code{TRUE} print the constructed
  command line instead of passing it to
  \code{\link{system}}.}
}
\description{
  Wrapper for the new NCBI BLAST+ tools
}
\keyword{internal}

