\name{blastReport-class}
\alias{blastReport-class}
\alias{show,blastReport-method}
\title{blastReport class}
\arguments{
  \item{...}{Slots for \sQuote{blastReport} instances.}
}
\description{
  blastReport is an S4 class that provides a container for
  data retrieved by calls to the NCBI Blast utility.
}
\details{
  blastReport objects have eleven slots: \describe{
  \item{program}{} \item{version}{} \item{reference}{}
  \item{db}{} \item{iter_num}{} \item{query}{}
  \item{hits}{} \item{params}{} \item{stats}{}
  \item{message}{} \item{data}{} }
}

