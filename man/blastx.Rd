\name{blastx}
\alias{blastx}
\title{Wrapper for the NCBI Translated Query-Protein Subject BLAST}
\usage{
  blastx(query, query_gencode=1, db="nr", out=NULL,
    outfmt="xml", max_hits=20, evalue=10,
    matrix="BLOSUM62", remote=FALSE, ...)
}
\arguments{
  \item{query}{Query sequences as path to a FASTA file, an
  \code{\linkS4class{XStringSet}} object, or a character
  vector.}

  \item{query_gencode}{Genetic code used to translate the
  query (default: 1).}

  \item{db}{The database to BLAST against (default: nr).}

  \item{out}{(optional) Output file for alignment. If
  \code{NULL} and the BLAST result is returned as a
  \code{\linkS4class{blastReport}} or
  \code{\linkS4class{blastTable}} object.}

  \item{outfmt}{Output format, \code{'xml'} or
  \code{'table'}.}

  \item{max_hits}{How many hits to return (default: 20).}

  \item{evalue}{Expect value cutoff (default: 10).}

  \item{matrix}{Scoring matrix name (default: BLOSUM62).}

  \item{remote}{Execute search remotely.}

  \item{...}{Additional parameters passed on to the BLAST
  commmand line tools. See
  \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
  for a description of common options.}
}
\description{
  Run \code{blastx()} without arguments to print usage and
  arguments description.
}
\examples{
##
}
\seealso{
  Other blast.applications: \code{\link{blastn}},
  \code{\link{blastn_short}}, \code{\link{blastp}},
  \code{\link{blastp_short}}, \code{\link{dc_megablast}},
  \code{\link{makeblasttdb}}, \code{\link{megablast}},
  \code{\link{qblast}}, \code{\link{rpsblast}},
  \code{\link{tblastn}}, \code{\link{tblastx}},
  \code{\link{update_blastdb}}
}

