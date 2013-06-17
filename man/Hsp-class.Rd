\name{Hsp-class}
\alias{Hsp-class}
\alias{HspList-class}
\title{Contains information parsed from an XML NCBI BLAST hsp element.}
\description{
  Information about each high-scoring pair (hsp) within a
  hit contained in a single \sQuote{Hsp} element.
}
\details{
  \section{Slots}{ \describe{ \item{\code{hsp_num}:}{The
  number of the hsp; \code{"integer"}.}
  \item{\code{score}:}{The BLAST score of the hsp;
  \code{"numeric"}.} \item{\code{bit_score}:}{The bit score
  of the hsp; \code{"numeric"}.} \item{\code{evalue}:}{The
  expect value; \code{"numeric"}.}
  \item{\code{identity}:}{Number of identities;
  \code{"integer"}.} \item{\code{positive}:}{Number of
  positives; \code{"integer"}.} \item{\code{gaps}:}{Number
  of gaps; \code{"integer"}.}
  \item{\code{query_from}:}{Start residue for query
  sequence; \code{"integer"}.} \item{\code{query_to}:}{End
  residue for query sequence; \code{"integer"}.}
  \item{\code{hit_from}:}{Start residue for hit sequence;
  \code{"integer"}.} \item{\code{hit_to}:}{End residue for
  hit sequence; \code{"integer"}.}
  \item{\code{query_frame}:}{\code{"integer"}.}
  \item{\code{hit_frame}:}{\code{"integer"}.}
  \item{\code{qseq}:}{Query sequence; \code{"XString"}.}
  \item{\code{hseq}:}{Hit sequence; \code{"XString"}.}
  \item{\code{match}:}{Match sequence/midline;
  \code{"XString"}.} \item{\code{query_env}:}{Shared
  container for \code{query_id}, \code{query_def},
  \code{query_len}, and \code{hit_len};
  \code{"environment"}.} } }
}

