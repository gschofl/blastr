\name{blastReport-class}
\alias{blastReport-class}
\alias{hits,blastReport-method}
\alias{show,blastReport-method}
\title{blastReport class}
\usage{
  .blastReport(...)
}
\arguments{
  \item{...}{Slots for \sQuote{blastReport} objects.}
}
\description{
  blastReport is an S4 class that provides a container for
  data retrived by calls to the NCBI Blast utility.
}
\details{
  blastReport objects have eight slots: \describe{
  \item{program}{} \item{version}{} \item{reference}{}
  \item{db}{} \item{query}{} \item{params}{} \item{stats}{}
  \item{hits}{} }
}

