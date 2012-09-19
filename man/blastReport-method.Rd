\name{query}
\alias{getAccn}
\alias{getAccn-methods}
\alias{getBitscore}
\alias{getBitscore-methods}
\alias{getDesc}
\alias{getDesc-methods}
\alias{getEvalue}
\alias{getEvalue-methods}
\alias{getHitFrame}
\alias{getHitFrame-methods}
\alias{getHitRange}
\alias{getHitRange-methods}
\alias{getHitSeq}
\alias{getHitSeq-methods}
\alias{getId}
\alias{getId-methods}
\alias{getQueryRange}
\alias{getQueryRange-methods}
\alias{getQuerySeq}
\alias{getQuerySeq-methods}
\alias{hits}
\alias{hits-methods}
\alias{index,blastReport-method}
\alias{query}
\alias{query-methods}
\title{Getter methods for blast records}
\usage{
  \S4method{index}{blastReport}(x, attributes = FALSE, ...)

  getId(x, db="gi")

  getDesc(x)

  getAccn(x)

  getBitscore(x)

  getEvalue(x)

  getQueryRange(x)

  getHitRange(x)

  getHitFrame(x)

  getQuerySeq(x)

  getHitSeq(x)
}
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

