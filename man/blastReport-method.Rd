\docType{methods}
\name{query}
\alias{getAccn}
\alias{getBitscore}
\alias{getDesc}
\alias{getEvalue}
\alias{getHitFrame}
\alias{getHitRange}
\alias{getHitSeq}
\alias{getId}
\alias{getQueryRange}
\alias{getQuerySeq}
\alias{hits}
\alias{index}
\alias{query}
\title{Getter methods for blast records}
\arguments{
  \item{x}{A \code{\linkS4class{blastRecord}} object.}

  \item{...}{Additional arguments}

  \item{x}{A \code{\linkS4class{blastRecord}} object.}

  \item{attributes}{Include definition as a name
  attribute.}

  \item{x}{A \code{\linkS4class{blastRecord}} object.}

  \item{...}{Additional arguments}

  \item{x}{A \code{\linkS4class{blastRecord}} or
  \code{\linkS4class{hit}} object.}

  \item{db}{Database tag (e.g.: 'gi', 'gb', 'emb', 'ref',
  'gnl', ...)}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}

  \item{x}{A \code{\link{blastReport-class}} or
  \code{\link{hit-class}} object.}
}
\value{
  Accession numbers as specified by \code{db} for each hit.

  Description for each hit.

  Accession number for each hit.

  Bitscores for each hit.

  Evalue for each hit.

  Start and end position for each query.

  Start and end position for each hit.

  Frame for each hit (1 or -1)

  Query sequence(s) as a
  \code{\link[Biostrings]{DNAStringSet}}

  Hit sequence(s) as a
  \code{\link[Biostrings]{DNAStringSet}}
}
\description{
  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records
}

