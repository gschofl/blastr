\docType{methods}
\name{query}
\alias{getAccn}
\alias{getAlignLen}
\alias{getBitscore}
\alias{getDesc}
\alias{getEvalue}
\alias{getGaps}
\alias{getGeneID}
\alias{getHitCov}
\alias{getHitDef}
\alias{getHitEnd}
\alias{getHitFrame}
\alias{getHitID}
\alias{getHitLen}
\alias{getHitNum}
\alias{getHitRange}
\alias{getHitSeq}
\alias{getHitStart}
\alias{getHSeq}
\alias{getHspHitID}
\alias{getHspID}
\alias{getHspNum}
\alias{getId}
\alias{getIdentity}
\alias{getMidline}
\alias{getPercIdent}
\alias{getPositive}
\alias{getQSeq}
\alias{getQueryCov}
\alias{getQueryDef}
\alias{getQueryEnd}
\alias{getQueryFrame}
\alias{getQueryLen}
\alias{getQueryRange}
\alias{getQuerySeq}
\alias{getQueryStart}
\alias{getScore}
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

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}

  \item{con}{A \code{\link{blastReportDB-class}} connection
  to a database}

  \item{x}{A query ID}
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

  Query definition(s) as a list

  Query length(s) as a list

  Hit id(s) as a list

  hit num(s) as a list

  Gene ID(s) of a hit as a list

  accession number(s) of a hit as a list

  hit definition(s) as a list

  hit length(s) as a list

  hit id(s) of a hsp as a list

  hsp id(s) as a list

  hsp num(s) as a list

  bitscore(s) as a list

  score(s) as a list

  evalue(s) as a list

  start position(s) of a hsp in query sequence as a list

  ending position(s) of a hsp in query sequence as a list

  start position(s) of a hsp in hit sequence as a list

  ending position(s) of a hsp in hit sequence as a list

  query frame(s) as a list

  hit frame(s) as a list

  identity(s) of a hsp as a list

  sum of positive(s) of a hsp as a list

  number off gap(s) as a list

  Length of Alignment(s) as a list

  query sequence(s) as a list

  hit sequence(s) as a list

  midline(s) as a list

  percent idetity(s) of a hsp as a list

  query coverage(s) of a hsp as a list

  hit coverage(s) of a hsp as a list
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

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records

  Getter methods for blast records
}

