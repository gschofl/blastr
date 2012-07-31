\name{qblast}
\alias{qblast}
\title{Do a BLAST search using the QBLAST URL API}
\usage{
  qblast(query,
    program = c("megablast", "blastn", "blastp", "rpsblast", "blastx", "tblastn", "tblastx"),
    db = "nr", outfmt = c("xml", "tabular", "html"),
    max_hits = 20, entrez_query = "", evalue = 10, ...,
    .params = list(), display = FALSE, update_time = 4)
}
\arguments{
  \item{query}{An accession number, a gi, or a fasta
  sequence provided either as an
  \code{\link[Biostrings]{XString}} or
  \code{\link[Biostrings]{XStringSet}} object, or as a
  named character vector.}

  \item{program}{The blast application to use. One of
  'megablast' (default), 'blastn', 'blastp', 'rpsblast',
  blastx', 'tblastn', 'tblastx'.}

  \item{db}{Blast database name. Defaults to 'nr'.}

  \item{outfmt}{Output format. One of 'xml', 'tabular', or
  'html'. See Details.}

  \item{max_hits}{The maximum number of hits to return.}

  \item{evalue}{Expectation value threshold for saving
  hits.}

  \item{entrez_query}{An Entrez query used to limit the
  results.}

  \item{...}{The name-value pairs of parameters passed on
  to the QBLAST URL API.}

  \item{.params}{A named list of parameters passed on to
  the QBLAST URL API.}

  \item{display}{Display the query result in a web
  browser.}

  \item{update_time}{How often to poll the blast server for
  results.}
}
\value{
  A \code{\link{blastReport-class}} object.
}
\description{
  Do a BLAST search using the QBLAST URL API
}
\examples{
##
}
\seealso{
  Other blast.applications: \code{\link{blastn}},
  \code{\link{blastn_short}}, \code{\link{blastp}},
  \code{\link{blastp_short}}, \code{\link{blastx}},
  \code{\link{dc_megablast}}, \code{\link{makeblasttdb}},
  \code{\link{megablast}}, \code{\link{tblastn}},
  \code{\link{tblastx}}
}

