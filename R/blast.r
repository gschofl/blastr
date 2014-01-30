#' @include blastReport-class.r
#' @importFrom RCurl getForm
#' @importFrom XML htmlParse
#' @importFrom Biostrings toString
#' @importFrom stringr str_match
NULL


#' Wrapper for NCBI makeblastdb
#' 
#' @param input_file Input file/database name. Multiple file/database names
#' can be provided as a character vector.
#' @param input_type Type of data specified in input file. One of \dQuote{fasta},
#' \dQuote{blastdb}, \dQuote{asn1_bin}, or \dQuote{asn1_txt}.
#' @param dbtype Molecule type of target db. (\dQuote{nucl} or \dQuote{prot})
#' @param ... further arguments passed to makeblastdb.
#' @param show_log print log file.
#' @param show_cmd print the command line instead of executing it.
#' 
#' @family blast applications
#' @export
makeblasttdb <- function(input_file, input_type = 'fasta', dbtype = 'nucl',
                         ..., show_log=TRUE, show_cmd=FALSE) {
  assert_that(has_command('makeblastdb'))
  if (missing(input_file)) {
    return(SysCall("makeblastdb", help=TRUE, redirection=FALSE))
  }
  ## assert that multiple input files are present and readable
  lapply(input_file, Compose(assert_that, is.readable))
  
  if (length(input_file) > 1) {
    input_file <- sprintf("\"%s\"", paste(input_file, collapse=" "))
  }
  
  input_type <- match.arg(input_type, c("fasta","blastdb","asn1_bin","asn1_txt"))
  dbtype <- match.arg(dbtype, c("nucl","prot"))
  
  o <- dots(...)
  if (!is.null(o$logfile)) {
    logfile <- o$logfile
  } else {
    logfile <- replace_ext(input_file[[1]], "log")
  }
  SysCall(exec="makeblastdb", infile=NULL, outfile=NULL,
          `in`=input_file, input_type=input_type,
          dbtype=dbtype, logfile=logfile, ..., style="unix",
          show_cmd=show_cmd)  
  if (show_log && assert_that(is.readable(logfile))) {
    cat(paste(readLines(logfile), collapse="\n"))
  }
}

#' Wrapper for update_blastdb.pl
#' 
#' Download pre-formatted BLAST databases from NCBI ftp site.
#' 
#' @param ... Blast db to download
#' @param destdir Destination directory for databases.
#' @param decompress if \code{TRUE}, decompresses the archives in \code{destdir}
#' and deletes the downloaded archives.
#' @param showall if \code{TRUE}, show all available pre-formatted BLAST
#' databases
#' @param passive Use passive FTP
#' @param timeout Timeout on connection to NCBI (default: 120 seconds)
#' @param force Force download even if there is a archive already on local
#' directory
#' 
#' @family blast applications
#' @export
update_blastdb <- function(..., destdir=getOption("blastr.blastdb.path") %||% '.',
                           decompress=TRUE, showall=FALSE, passive=FALSE,
                           timeout=120, force=FALSE) {
  destdir <- normalizePath(destdir)
  assert_that(has_command("update_blastdb"))
  blastdb <- unlist(dots(...))
  args <- list(decompress=FALSE, passive=passive, force=force, timeout=timeout)  
  if (showall) {
    ans <- SysCall('update_blastdb', showall=TRUE, style='gnu', intern=TRUE)
    return( ans[-1] )
  } else if (all_empty(blastdb) && all(are_false(args)[-4])) {
    SysCall('update_blastdb')
  } else {
    available_dbs <- SysCall('update_blastdb', showall=TRUE,
                             style='gnu', intern=TRUE)[-1]
    if (!all(idx <- blastdb %in% available_dbs)) {
      stop(sQuote(paste0(blastdb[!idx], collapse=", ")),
           " not among the available databases")
    }
    cwd <- setwd(destdir)
    on.exit(setwd(cwd))
    SysCall('update_blastdb', verbose=TRUE, stdin=paste0(blastdb, collapse=' '),
            args=args, style="gnu", redirection=FALSE, show_cmd=FALSE,
            intern=FALSE)
  } 
  if (decompress) {
    assert_that(has_command("tar"))
    archives <- dir(destdir, pattern="gz$", full.names=TRUE)
    if (length(archives) > 0) {
      sapply(archives, function (a) {
        system(paste("tar zxvpf", a))
        unlink(a)
      })
    }
  }
}


#' Wrapper for the new NCBI BLAST+ tools
#' 
#' @param program One of blastn, blastp, blastx, tblastn, or tblastx
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db Blast database name (defaults to 'nt' for blastn and 'nr' for
#' the rest).
#' @param outfmt XML or table (default: XML)
#' @param max_hits Maximum number of hits  to return. (defaults: 20)
#' Sets the '-max_target_seqs' parameter internally.
#' @param strand Query strand(s) to seach against database.
#' @param ... Arguments passed on to the blast commmand line tools.
#' @param intern Set \code{TRUE} if no '-out' argument is specified.
#' Captures the blast output in an R character vector.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' @param parse
#' 
#' @keywords internal
.blast <- function(exec, query, db, outfmt = 'xml', max_hits = 20,
                   strand = 'both', ..., intern = FALSE, show_cmd = FALSE,
                   parse = TRUE) {
  blastdb.old <- Sys.getenv("BLASTDB")
  on.exit(Sys.setenv(BLASTDB = blastdb.old))
  exec <- match.arg(exec, c("blastn", "blastp", "blastx", "tblastn",
                            "tblastx", "rpsblast+", 'rpstblastn'))
  assert_that(has_command(exec))
  strand <- match.arg(strand, c("both", "plus", "minus"))
  outfmt <- switch(match.arg(outfmt, c('xml', 'table')), xml=5, table=7)
  dot_args <- list(...)
  
  # dealing with query
  if (missing(query)) {
    return(SysCall(exec, help=TRUE, intern=FALSE))
  }
  if (exec %in% c("blastp", "blastp_short", "rpsblast+")) {
    strand <- NULL
  }
  
  inp <- make_blast_query(query)
  query <- inp[["query"]]
  on.exit(unlink(query))
  nqueries <- inp[["nqueries"]]
  parse_deflines <- inp[["parse_deflines"]]
  
  # set a number of defaults different from the internal defaults of
  # the blast applications
  if (missing(db)) {
    db <- switch(exec, blastn='nt', `rpsblast+`='Cdd', 'nr')
  }
  
  ## if we query the NCBI blast server we don't want to fetch the path to local
  ## blast dbs
  if (!(dot_args$remote %||% FALSE)) {
    ## check if the supplied path to a BLAST db is valid or if a valid path is
    ## set by the global option blastr.blastdb.path
    db <- set_blastdb(db)
  }
  args <- merge_list(
    dot_args,
    list(query=query, db=db, outfmt=outfmt,
         num_descriptions=NULL, num_alignments=NULL,
         max_target_seqs=max_hits, strand=strand,
         parse_deflines=parse_deflines)
  ) 
  # check if 'out' was specified, otherwise return results internally
  intern <- if (is.null(args[["out"]])) TRUE else FALSE
  cat(paste0("Blasting ", nqueries, " queries [", exec, "]", "\n"), sep="")
  res <- SysCall(exec=exec, args=args, redirection= FALSE, style="unix",
                 show_cmd=show_cmd, intern=intern)
  if (intern && parse && !has_attr(res, 'status')) {
    if (outfmt == 5) {
      blastReport(res, asText=TRUE)
    } else if (outfmt == 7) {
      blastTable(res)
    } else {
      res
    }
  } else res 
}


set_blastdb <- function(db) {
  db.ori <- db
  if (!is.blastdb(db)) {
    db <- normalizePath(file.path(getOption("blastr.blastdb.path"), db), mustWork = FALSE)
    if (!is.blastdb(db)) {
      stop("'", basename(db.ori), "' is not a valid BLAST database", call. = FALSE)
    }
  }
  Sys.setenv(BLASTDB = dirname(db))
  basename(db)
}


is.blastdb <- function(db) {
  ## Protein and DNA database estensions, respectively.
  ## three of these must be present
  db_ext <- c("pin", "psq", ".phr", "nin", "nsq", "nhr")
  if (length(dbs <- dir(dirname(db), pattern = paste0(basename(db), '\\..+'))) > 0) {
    sum(db_ext %in% strsplitN(dbs, "\\.", 1, "end")) == 3
  } else {
    FALSE
  }
}


#' Wrapper for the NCBI Nucleotide-Nucleotide BLAST
#' 
#' \itemize{
#' \item{\code{blastn}} is the traditional BLASTN requiring an exact
#' match of 11.
#' \item{\code{blastn_short}} is BLASTN optimised for sequences shorter
#' than 50 bases.
#' \item{\code{megablast}} is the traditional megablast used to find
#' very similar sequences.
#' \item{\code{dc_megablast}} is discontiguous megablast used to find
#' more distant (e.g. interspecies) sequences.
#' }
#' 
#' Run \code{blastn()} without arguments to print usage and
#' arguments description.
#' 
#' @usage blastn(query, db="nt", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, remote=FALSE, ...)
#'
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db The database to BLAST against (default: nt).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export blastn blastn_short megablast dc_megablast
#' @aliases blastn blastn_short megablast dc_megablast
#' @examples
#' ##
blastn <- Partial(.blast, exec = "blastn", task = "blastn")

#' @usage blastn_short(query, db="nt", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, remote=FALSE, ...)
#' @export
#' @rdname blastn
blastn_short <- Partial(.blast, exec = "blastn", task = "blastn-short")

#' @usage megablast(query, db="nt", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, remote=FALSE, ...)
#' @export
#' @rdname blastn
megablast <- Partial(.blast, exec = "blastn", task = "megablast")

#' @usage dc_megablast(query, db="nt", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, remote=FALSE, ...)
#' @export
#' @rdname blastn
dc_megablast <- Partial(.blast, exec = "blastn", task = "dc-megablast")

#' Wrapper for the NCBI Protein-Protein BLAST
#' 
#' \code{blastp} is the traditional BLASTP.
#' \code{blastp_short} is BLASTP optimised for residues shorter than 30.
#' 
#' Run \code{blastp()} without arguments to print usage and
#' arguments description.
#' 
#' @usage blastp(query, db="nr", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, matrix="BLOSUM62", remote=FALSE, ...)
#'
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db The database to BLAST against (default: nr).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param matrix Scoring matrix name (default: BLOSUM62).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export blastp blastp_short
#' @aliases blastp blastp_short
#' @examples
#' ##
blastp <- Partial(.blast, exec = "blastp", task = "blastp")

#' @usage blastp_short(query, db="nr", out=NULL, outfmt="xml", max_hits=20,
#'  evalue=10, matrix="BLOSUM62", remote=FALSE, ...)
#' @export
#' @rdname blastp
#' @inheritParams blastp
blastp_short <- Partial(.blast, exec = "blastp", task = "blastp-short")

#' Wrapper for the NCBI Translated Query-Protein Subject BLAST
#' 
#' Run \code{blastx()} without arguments to print usage and
#' arguments description.
#' 
#' @usage blastx(query, query_gencode=1, db="nr", out=NULL, outfmt="xml",
#'  max_hits=20, evalue=10, matrix="BLOSUM62", remote=FALSE, ...)
#' 
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param query_gencode Genetic code used to translate the query (default: 1).
#' @param db The database to BLAST against (default: nr).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param matrix Scoring matrix name (default: BLOSUM62).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export blastx
#' @aliases blastx
#' @examples
#' ##
blastx <- Partial(.blast, exec = "blastx")

#' Wrapper for the NCBI Translated Query-Protein Subject BLAST
#' 
#' Run \code{tblastx()} without arguments to print usage and
#' arguments description.
#' 
#' @usage tblastx(query, query_gencode=1, db="nr", out=NULL, outfmt="xml",
#'  max_hits=20, evalue=10, matrix="BLOSUM62", remote=FALSE, ...)
#' 
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param query_gencode Genetic code used to translate the query (default: 1).
#' @param db The database to BLAST against (default: nr).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param matrix Scoring matrix name (default: BLOSUM62).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export tblastx
#' @aliases tblastx
#' @examples
#' ##
tblastx <- Partial(.blast, exec = "tblastx")

#' Wrapper for the NCBI Protein Query-Translated Subject BLAST
#' 
#' Run \code{tblastn()} without arguments to print usage and
#' arguments description.
#' 
#' @usage tblastn(query, db="nr", out=NULL, outfmt="xml", max_hits=20,
#' evalue=10, matrix="BLOSUM62", remote=FALSE, ...)
#' 
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db The database to BLAST against (default: nr).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param matrix Scoring matrix name (default: BLOSUM62).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export tblastn
#' @aliases tblastn
#' @examples
#' ##
tblastn <- Partial(.blast, exec = "tblastn")


#' Wrapper for the NCBI Reversed Position Specific Blast
#' 
#' Run \code{rpsblast()} without arguments to print usage and
#' arguments description.
#' 
#' @usage rpsblast(query, db="Cdd", out=NULL, outfmt="xml", max_hits=20,
#' evalue=10, remote=FALSE, ...)
#' 
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db The database to BLAST against (default: Cdd).
#' @param out (optional) Output file for alignment.
#' If \code{NULL} and the BLAST result is returned as
#' a \code{\linkS4class{blastReport}} or \code{\linkS4class{blastTable}}
#' object.
#' @param outfmt Output format, \code{'xml'} or \code{'table'}.
#' @param max_hits How many hits to return (default: 20).
#' @param evalue Expect value cutoff (default: 10).
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' 
#' @family blast applications
#' @export rpsblast
#' @aliases rpsblast
#' @examples
#' ##
rpsblast <- Partial(.blast, exec = "rpsblast+")


#' Do a BLAST search using the QBLAST URL API
#' 
#' @param query Query sequences as path to a FASTA file,
#' an \code{\linkS4class{XStringSet}} object, or a character vector.
#' @param db The database to BLAST against.
#' @param program The blast application to use. One of 'megablast'
#' (default), 'blastn', 'blastp', 'rpsblast', blastx', 'tblastn', 'tblastx'.
#' @param db Blast database name. (defaults: 'nt').
#' @param outfmt Output format. One of 'xml' or 'tabular'.
#' See Details.
#' @param max_hits The maximum number of hits to return.
#' @param evalue Expectation value threshold for saving hits.
#' @param entrez_query An Entrez query used to limit the results. 
#' @param ... The name-value pairs of parameters passed on to the QBLAST
#' URL API.
#' @param .params A named list of parameters passed on to the QBLAST
#' URL API.
#' @param update_time How often to poll the blast server for results.
#' 
#' @return A \code{\link{blastReport-class}} object.
#'
#' @family blast applications
#' @export
#' @examples
#' ##
qblast <- function(query, exec = 'megablast', db = 'nt', outfmt = 'xml',
                   max_hits = 20, entrez_query = '', evalue = 10, ...,
                   .params = list(), update_time = 4)
{
  if (missing(query))
    stop("No query provided")
  
  qbd <- megablast <- rpsblast <- FALSE
  exec <- match.arg(exec, c("megablast","blastn","blastp", "rpsblast",
                                  "blastx","tblastn", "tblastx"))
  outfmt <- match.arg(outfmt, c("xml", "tabular"))
  
  if (exec == "megablast") {
    exec <- "blastn"
    megablast <- TRUE
  }
  
  if (exec == "rpsblast") {
    exec <- "blastp"
    rpsblast <- TRUE
  }
  
  if (outfmt == "xml") {
    format_type <- "XML"
    alignment_view <- "Pairwise"
  } 
  else if (outfmt == "tabular") {
    format_type <- "Text"
    alignment_view <- "Tabular"
  }
  
  if (is(query, "XStringSet") || is(query, "XString")) {
    query <- paste0(">", names(query)[[1L]], "\n", Biostrings::toString(query[[1L]]))
    qbd <- TRUE
  } else if (!is.null(names(query))) {
    query <- paste0(">", names(query)[[1L]], "\n", as.character(query[[1L]]))
    qbd <- TRUE
  } else {
    query <- as.character(query)
  }
  
  if (nchar(entrez_query) > 0L) {
    .params <- merge_list(.params, list(entrez_query=.escape(entrez_query)))
  }
  
  .params <- merge_list(list(query=query,
                             program=exec,
                             database=db,
                             expect=evalue,
                             hitlist_size=max_hits,
                             query_believe_defline=if (qbd) "true" else "false",
                             megablast=if (megablast) "yes" else "no",
                             service=if (rpsblast) "rpsblast" else "plain",
                             client="web",
                             cmd="Put"),
                        list(...),
                        .params)
  
  names(.params) <- toupper(names(.params)) 
  base_url <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"   
  response <- getForm(base_url, .params=.params, style="post")
  # parse out the request id and the estimated time to completion
  r <- xvalue(htmlParse(response),
              "/descendant::comment()[contains(.,'RID =')]")
  rid <- str_match(r, "RID = (.[^(\n)]*)")[,2]
  
  # poll for results
  while (TRUE) {
    
    Sys.sleep(update_time)
    cat(sprintf(" Searching ... (update in %s seconds)\n", update_time))
    response <- getForm(uri=base_url,
                        RID=rid,
                        FORMAT_OBJECT="SearchInfo",
                        CMD="Get",
                        style="post") 
    status <- xvalue(htmlParse(response),
                     "/descendant::comment()[contains(.,'Status=')]")
    status <- str_match(status, "Status=(.[^(\n)]*)")[,2]
    cat(sprintf("Status: %s\n", status))
    
    if (status == "WAITING") {
      next
    }
    
    if (status == "FAILED") {
      stop(sprintf("Search %s failed\n", rid))
    }
    
    if (status == "UNKNOWN") {
      stop(sprintf("Search %s expired\n", rid))
    }
    
    if (status == "READY") {
      there_are_hits <- xvalue(htmlParse(response),
                               "/descendant::comment()[contains(.,'ThereAreHits=')]")
      there_are_hits <- str_match(there_are_hits, "ThereAreHits=(.[^(\n)]*)")[,2]
      if (there_are_hits == "yes") {
        break
      }
      else {
        cat("No hits found\n")
        break
      }
    }
  }
  
  # retrieve results
  response <- getForm(uri=base_url,
                      RID=rid,
                      FORMAT_TYPE=format_type,
                      FORMAT_OBJECT="Alignment",
                      DESCRIPTIONS=max_hits,
                      ALIGNMENTS=max_hits,
                      ALIGNMENT_VIEW=alignment_view,
                      NEW_VIEW="true",
                      DISPLAY_SORT="2",
                      NCBI_GI="true",
                      CMD="Get",
                      style="post")
  
  if (outfmt == "tabular") {
    response <- xvalue(htmlParse(response), "//pre")
  }
  
  response
}
