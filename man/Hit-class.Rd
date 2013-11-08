\name{Hit-class}
\alias{Hit-class}
\alias{HitList-class}
\title{Hit-class}
\usage{
  new_Hit(...)
}
\description{
  A container for information parsed from an XML NCBI BLAST
  hit element. Information about multiple Blast hits within
  a query is contained in a \code{HitList}.
}
\details{
  \describe{ \item{\code{hit_num}:}{The number of the hit;
  \code{"integer"}.} \item{\code{hit_def}:}{Hit definition;
  \code{"DeflineSet"}.} \item{\code{hit_acc}:}{Accession
  number; \code{"character"}.}
  \item{\code{hit_len}:}{Length of hit; \code{"integer"}.}
  \item{\code{hsps}:}{List of HSPs;
  \code{"\linkS4class{HspList}"}.}
  \item{\code{query_env}:}{Shared container for
  \code{query_id}, \code{query_def}, \code{query_len}, and
  \code{hit_len}; \code{"environment"}.} }
}

