\name{makeblasttdb}
\alias{makeblasttdb}
\title{Wrapper for NCBI makeblastdb}
\usage{
  makeblasttdb(input_file,
    input_type = c("fasta", "blastdb", "asn1_bin", "asn1_txt"),
    dbtype = c("nucl", "prot"), ..., show_log = TRUE,
    show_cmd = FALSE)
}
\arguments{
  \item{input_file}{Input file/database name.}

  \item{input_type}{Type of data specified in input file.}

  \item{dbtype}{Molecule type of target db.}

  \item{...}{further arguments passed to makeblastdb.}

  \item{show_log}{print log file.}

  \item{show_cmd}{print the command line instead of
  executing it.}
}
\description{
  Wrapper for NCBI makeblastdb
}
\seealso{
  Other blast.applications: \code{\link{blastn}},
  \code{\link{blastn_short}}, \code{\link{blastp}},
  \code{\link{blastp_short}}, \code{\link{blastx}},
  \code{\link{dc_megablast}}, \code{\link{megablast}},
  \code{\link{qblast}}, \code{\link{tblastn}},
  \code{\link{tblastx}}
}

