\name{BlastParameters}
\alias{BlastParameters}
\alias{BlastParameters-class}
\title{blastParameters-class}
\description{
  A container for blast parameters and statistics.
}
\details{
  \describe{ \item{\code{program}:}{The BLAST flavour that
  generated the data; \code{"character"}.}
  \item{\code{matrix}:}{Name of the matrix (\code{NA} for
  nucleotide blast); \code{"character"}.}
  \item{\code{expect}:}{Cutoff value; \code{"numeric"}.}
  \item{\code{penalties}:}{Open and extend penalties;
  \code{"numeric"}.} \item{\code{sc_match}:}{Match score
  for nucleotide-nucleotide comparison; \code{"integer"}.}
  \item{\code{sc_mismatch}:}{Mismatch penalty for
  nucleotide-nucleotide comparison; \code{"integer"}.}
  \item{\code{filter}:}{Filter string; \code{"character"}.}
  \item{\code{num_sequences}:}{Number of sequences in the
  database; \code{"character"}.}
  \item{\code{num_letters}:}{Number of letters in the
  database; \code{"character"}.}
  \item{\code{hsp_length}:}{Effective HSP length;
  \code{"numeric"}.}
  \item{\code{effective_space}:}{Effective search space;
  \code{"numeric"}.} \item{\code{ka_params}:}{kappa,
  lambda, entropy; \code{"numeric"}.} }
}
\keyword{internal}

