\name{blastReportDB-class}
\alias{blastReportDB-class}
\title{blastReportDB-class}
\usage{
  .blastReportDB(...)
}
\description{
  \sQuote{\bold{blastReportDB}}: A connection to an SQLite
  database containing blast records organised in three
  tables:
}
\details{
  \bold{query} with fields:

  \itemize{ \item query_id INTEGER Primary key \item
  query_def TEXT \item query_len INTEGER }

  \bold{hit} with fields:

  \itemize{ \item query_id INTEGER \item hit_id INTEGER
  Primary key \item hit_num INTEGER \item gene_id TEXT
  \item accession TEXT \item definition TEXT \item length
  INTEGER }

  \bold{hsp} with fields:

  \itemize{ \item query_id INTEGER \item hit_id INTEGER
  \item hsp_id INTEGER Primary key \item hsp_num INTEGER
  \item bit_score FLOAT \item score INTEGER \item evalue
  FLOAT \item query_from INTEGER \item query_to INTEGER
  \item hit_from INTEGER \item hit_to INTEGER \item
  query_frame INTEGER \item query_frame INTEGER \item
  identity INTEGER \item positive INTEGER \item gaps
  INTEGER \item align_len INTEGER \item qseq TEXT \item
  hseq TEXT \item midline TEXT }
}
\seealso{
  The constructors \code{\link{blastReportDB}} and
  \code{\link{blastReportDBConnect}}, and the BLAST classes
  \code{\linkS4class{blastReport}} and
  \code{\linkS4class{blastTable}}.
}

