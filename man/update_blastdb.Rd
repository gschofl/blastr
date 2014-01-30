\name{update_blastdb}
\alias{update_blastdb}
\title{Wrapper for update_blastdb.pl}
\usage{
update_blastdb(..., destdir = getOption("blastr.blastdb.path") \%||\% ".",
  decompress = TRUE, showall = FALSE, passive = FALSE, timeout = 120,
  force = FALSE)
}
\arguments{
  \item{...}{Blast db to download}

  \item{destdir}{Destination directory for databases.}

  \item{decompress}{if \code{TRUE}, decompresses the
  archives in \code{destdir} and deletes the downloaded
  archives.}

  \item{showall}{if \code{TRUE}, show all available
  pre-formatted BLAST databases}

  \item{passive}{Use passive FTP}

  \item{timeout}{Timeout on connection to NCBI (default:
  120 seconds)}

  \item{force}{Force download even if there is a archive
  already on local directory}
}
\description{
Download pre-formatted BLAST databases from NCBI ftp site.
}
\seealso{
Other blast.applications: \code{\link{blastn}},
\code{\link{blastn}}, \code{\link{blastn_short}},
\code{\link{blastn_short}}, \code{\link{dc_megablast}},
\code{\link{dc_megablast}}, \code{\link{megablast}},
\code{\link{megablast}}; \code{\link{blastp}},
\code{\link{blastp}}, \code{\link{blastp_short}},
\code{\link{blastp_short}}; \code{\link{blastx}},
\code{\link{blastx}}; \code{\link{makeblasttdb}};
\code{\link{qblast}}; \code{\link{rpsblast}},
\code{\link{rpsblast}}; \code{\link{tblastn}},
\code{\link{tblastn}}; \code{\link{tblastx}},
\code{\link{tblastx}}
}

