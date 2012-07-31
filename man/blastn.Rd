\name{blastn}
\alias{blastn}
\alias{blastn_short}
\alias{dc_megablast}
\alias{megablast}
\title{Wrapper for the NCBI Nucleotide-Nucleotide BLAST}
\usage{
  blastn(query, db="nt", out=NULL, outfmt=5, max_hits=20,
    evalue=10, html=FALSE, remote=FALSE, ...,
    show_cmd=FALSE)

  blastn_short(query, db="nt", out=NULL, outfmt=5,
    max_hits=20, evalue=10, html=FALSE, remote=FALSE, ...,
    show_cmd=FALSE)

  megablast(query, db="nt", out=NULL, outfmt=5,
    max_hits=20, evalue=10, html=FALSE, remote=FALSE, ...,
    show_cmd=FALSE)

  dc_megablast(query, db="nt", out=NULL, outfmt=5,
    max_hits=20, evalue=10, html=FALSE, remote=FALSE, ...,
    show_cmd=FALSE)
}
\arguments{
  \item{query}{The sequence to search with provided as a
  file name, an \code{\link[Biostrings]{XString}} or
  \code{\link[Biostrings]{XStringSet}} object, or as a
  character vector.}

  \item{db}{The database to BLAST against.}

  \item{out}{Output file for alignment. If \code{NULL} the
  BLAST result is returned as an R character vector.}

  \item{outfmt}{Output format, Integer 1-11. Default is 5
  for XML output. Other options include 6 or 7 for tabular
  output and 0 for the traditional pairwise alignment
  view.}

  \item{max_hits}{How many hits to return.}

  \item{evalue}{Expectation value cutoff.}

  \item{html}{Produce HTML output.}

  \item{remote}{Execute search remotely.}

  \item{...}{Additional parameters passed on to the BLAST
  commmand line tools. See
  \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
  for a description of common options.}

  \item{show_cmd}{If \code{TRUE} print the constructed
  command line instead of passing it to
  \code{\link{system}}.}
}
\description{
  \itemize{ \item{\code{blastn}} is the traditional BLASTN
  requiring an exact match of 11.
  \item{\code{blastn_short}} is BLASTN optimised for
  sequences shorter than 50 bases. \item{\code{megablast}}
  is the traditional megablast used to find very similar
  sequences. \item{\code{dc_megablast}} is discontiguous
  megablast used to find more distant (e.g. interspecies)
  sequences. }
}
\details{
  Run \code{blastn()} without arguments to print usage and
  arguments description.
}
\examples{
##
}
\seealso{
  Other blast.applications: \code{\link{blastp}},
  \code{\link{blastp_short}}, \code{\link{blastx}},
  \code{\link{makeblasttdb}}, \code{\link{qblast}},
  \code{\link{tblastn}}, \code{\link{tblastx}}
}

