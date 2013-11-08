\name{Iteration-class}
\alias{Iteration-class}
\title{Iteration-class}
\description{
  An S4 class that holds data parsed from an NCBI BLAST XML
  \sQuote{Iteration} element. Each Iteration contains data
  about a query.
}
\section{Slots}{
  \describe{ \item{\code{iter_num}:}{The number of the
  iteration; \code{"integer"}.}
  \item{\code{query_id}:}{Query identifier;
  \code{"character"}.} \item{\code{query_def}:}{Query
  definition; \code{"character"}.}
  \item{\code{query_len}:}{Query length; \code{"integer"}.}
  \item{\code{hits}:}{List of hits;
  \code{"\linkS4class{HitList}"}.}
  \item{\code{query_env}:}{Shared container for
  \code{query_id}, \code{query_def}, \code{query_len}, and
  \code{hit_len}; \code{"environment"}.} }
}

