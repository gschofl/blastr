\name{blastTable-class}
\alias{blastTable-class}
\title{blastTable class}
\arguments{
  \item{...}{Slots for \sQuote{blastTable} instances.}
}
\description{
  blastTable is an S4 class that provides a container for
  data retrived by calls to the NCBI Blast utility.
}
\details{
  blastReport objects have ten slots: \describe{
  \item{program}{} \item{version}{} \item{reference}{}
  \item{db}{} \item{bit_score}{} \item{evalue}{}
  \item{mlog.evalue}{} \item{accession}{} \item{gi}{}
  \item{table}{} }
}

