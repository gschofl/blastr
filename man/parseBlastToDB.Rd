\name{parseBlastToDB}
\alias{parseBlastToDB}
\title{Parse large BLAST Reports into SQLite DB}
\usage{
  parseBlastToDB(blastFile, dbPath = "blast.db",
    max_hit = NULL, max_hsp = NULL, reset_at = 1000)
}
\arguments{
  \item{blastFile}{an XML BLAST Report.}

  \item{dbPath}{Path to database.}

  \item{max_hit}{How many hits should be parsed (default:
  all available)}

  \item{max_hsp}{How many hsps should be parsed (default:
  all available)}

  \item{reset_at}{After how many iterations should the
  parser dump the data into the db before continuing.}
}
\value{
  a connection object
}
\description{
  Parse large BLAST Reports into SQLite DB
}

