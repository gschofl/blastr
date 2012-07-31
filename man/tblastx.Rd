\name{tblastx}
\alias{tblastx}
\title{Wrapper for the NCBI Translated Query-Protein Subject BLAST}
\usage{
  tblastx(query, query_gencode=1, db="nr", out=NULL,
    outfmt=5, max_hits=20, evalue=10, matrix="BLOSUM62",
    html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
}
\arguments{
  \item{query}{The sequence to search with provided as a
  file name, an \code{\link[Biostrings]{XString}} or
  \code{\link[Biostrings]{XStringSet}} object, or as a
  character vector.}

  \item{query_gencode}{Genetic code used to translate the
  query (default: 1).}

  \item{db}{The database to BLAST against.}

  \item{out}{Output file for alignment. If \code{NULL} the
  BLAST result is returned as an R character vector.}

  \item{outfmt}{Output format, Integer 1-11. Default is 5
  for XML output. Other options include 6 or 7 for tabular
  output and 0 for the traditional pairwise alignment
  view.}

  \item{max_hits}{How many hits to return.}

  \item{evalue}{Expectation value cutoff.}

  \item{matrix}{Scoring matrix name (Default: BLOSUM62).}

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
  Run \code{tblastx()} without arguments to print usage and
  arguments description.
}
\examples{
##
}
\seealso{
  Other blast.applications: \code{\link{blastn}},
  \code{\link{blastn_short}}, \code{\link{blastp}},
  \code{\link{blastp_short}}, \code{\link{blastx}},
  \code{\link{dc_megablast}}, \code{\link{makeblasttdb}},
  \code{\link{megablast}}, \code{\link{qblast}},
  \code{\link{tblastn}}
}

