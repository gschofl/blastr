#' @include blastReport-class.r
#' @importFrom RCurl getForm
#' @importFrom XML htmlParse
#' @importFrom rmisc SysCall
#' @importFrom rmisc Curry
#' @importFrom rmisc merge_list
#' @importFrom rmisc has_command
#' @importFrom Biostrings toString
#' @importFrom stringr str_match
#' @importFrom assertthat has_attr
NULL


#' Wrapper for NCBI makeblastdb
#' 
#' @param input_file Input file/database name. Multiple file/database names
#' can be provided as a character vector.
#' @param input_type Type of data specified in input file.
#' @param dbtype Molecule type of target db. (\sQuote{nucl} or \sQuote{prot})
#' @param ... further arguments passed to makeblastdb.
#' @param show_log print log file.
#' @param show_cmd print the command line instead of executing it.
#' 
#' @family blast applications
#' @export
makeblasttdb <- function(input_file, input_type = 'fasta', dbtype = 'nucl',
                         ..., show_log=TRUE, show_cmd=FALSE) {
  assert_that(has_command('makeblastdb'))
  
  if (missing(input_file))
    return(SysCall("makeblastdb", help=TRUE, redirection=FALSE))

  ## assert that multiple input files are present
  lapply(input_file, assert_that %.% is.readable)
  
  if (length(input_file) > 1) {
    input_file <- sprintf("\"%s\"", paste(input_file, collapse=" "))
  }
  
  input_type <- match.arg(input_type, c("fasta","blastdb","asn1_bin","asn1_txt"))
  dbtype <- match.arg(dbtype, c("nucl","prot"))
  
  o <- list(...)
  if (!is.null(o$logfile))
    logfile <- o$logfile
  else
    logfile <- replace_ext(input_file[[1]], "log")

  SysCall(exec="makeblastdb", infile=NULL, outfile=NULL,
          `in`=input_file, input_type=input_type,
          dbtype=dbtype, logfile=logfile, ..., style="unix",
          show_cmd=show_cmd)
  
  if (show_log && assert_that(is.readable(logfile)))
    cat(paste(readLines(logfile), collapse="\n"))
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
#' @force force Force download even if there is a archive already on local
#' directory
#' 
#' @family blast applications
#' @export
update_blastdb <- function(..., destdir=".", decompress=FALSE, showall=FALSE,
                           passive=FALSE, timeout=120, force=FALSE) {
  assert_that(has_command("update_blastdb"))
  args <- list(decompress=decompress, passive=passive, force=force, timeout=timeout)
  if (showall) {
    ans <- SysCall('update_blastdb', showall=TRUE, style='gnu', intern=TRUE)
    return( ans[-1] )
  } else if (all(are_false(args)[-4])) {
    SysCall('update_blastdb')
  } else {
    available_dbs <- SysCall('update_blastdb', showall=TRUE,
                             style='gnu', intern=TRUE)[-1]
    blastdb <- c(...)
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
}


#' Wrapper for the new NCBI BLAST+ tools
#' 
#' @param program One of blastn, blastp, blastx, tblastn, or tblastx
#' @param query Query provided as an input file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param db Blast database name (defaults to 'nt' for blastn and 'nr' for
#' the rest).
#' @param outfmt Output format (defaults to XML output).
#' @param max_hits Maximum number of hits  to return. Defaults to 20.
#' Sets the '-max_target_seqs' parameter internally.
#' @param strand Query strand(s) to seach against database.
#' @param ... Arguments passed on to the blast commmand line tools.
#' @param intern Set \code{TRUE} if no '-out' argument is specified.
#' Captures the blast output in an R character vector.
#' @param input Used to pass a character vector to the standard input.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @keywords internal
blast <- function (program = 'blastp', query, db, outfmt, max_hits,
                   strand = 'both', ..., intern = FALSE, input = NULL,
                   show_cmd = FALSE, parse = TRUE) {
  
  program <- match.arg(program, c("blastn", "blastp", "blastx", "tblastn", "tblastx"))
  assert_that(has_command(program))
  strand <- match.arg(strand, c("both", "plus", "minus"))
  
  # dealing with query
  if (missing(query))
    return( SysCall(program, help=TRUE, intern=FALSE) )

  # set strand=NULL and translate sequence if blastp
  if (program %in% c("blastp","blastp_short")) {
    strand <- NULL
    transl <- TRUE
  } else {
    transl <- FALSE
  }
  
  inp <- make_blast_query(x=query, transl)
  input <- inp[["input"]]
  query <- inp[["query"]]
  deflines <- inp[["deflines"]]
  parse_deflines <- inp[["parse_deflines"]]

  # set a number of defaults different from the internal defaults of
  # the blast applications
  if (missing(db)) {
    if (program == "blastn")
      db <- "nt"
    else
      db <- "nr"
  }
  
  if (missing(outfmt))
    outfmt <- 5
  
  # always use tables with headers
  if (outfmt == 6)
    outfmt <- 7
  
  if (missing(max_hits))
    max_hits <- 20
  
  ## max_target_seqs is incompatible with output options 0 - 4
  if (as.integer(outfmt) <= 4) {
    num_descriptions <- max_hits
    num_alignments <- max_hits
    max_target_seqs <- NULL
  } else if (as.integer(outfmt) >= 5) {
    num_descriptions <- NULL
    num_alignments <- NULL
    max_target_seqs <- max_hits
  }
 
  args <- merge_list(list(...),
                     list(query=query,
                          db=db,
                          outfmt=outfmt,
                          num_descriptions=num_descriptions,
                          num_alignments=num_alignments,
                          max_target_seqs=max_target_seqs,
                          strand=strand,
                          parse_deflines=parse_deflines))

  # remove stdin', 'stdout' from the arguments list
  stdin <- args[["stdin"]]
  args[["stdin"]] <- NULL
  stdout <- args[["stdout"]]
  args[["stdout"]] <- NULL
  
  # check if 'out' is specified, otherwise return results internally
  intern <- if (is.null(args[["out"]])) TRUE else FALSE
  cat(paste0("Blasting ", length(deflines), " queries [", program, "]", "\n"), sep="")
  res <- SysCall(exec=program, args=args,
                 stdin=stdin, stdout=stdout,
                 redirection=if (is.null(stdin) && is.null(stdout)) FALSE else TRUE,
                 style="unix", show_cmd=show_cmd,
                 intern=intern, input=input)
  if (intern && outfmt == 5 && parse && !has_attr(res, 'status')) 
    return( blastReport(res, asText=TRUE) )
  else
    return( res )
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
#' @usage blastn(query, db="nt", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
#'
#' @param query The sequence to search with provided as a file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param db The database to BLAST against.
#' @param out Output file for alignment. If \code{NULL} and \code{outfmt=5}
#' (XML) the BLAST result is returned as \code{\linkS4class{blastReport}}.
#' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
#' Other options include 6 or 7 for tabular output and 0 for the
#' traditional pairwise alignment view.
#' @param max_hits How many hits to return.
#' @param evalue Expectation value cutoff.
#' @param html Produce HTML output.
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @family blast applications
#' @export blastn blastn_short megablast dc_megablast
#' @aliases blastn blastn_short megablast dc_megablast
#' @examples
#' ##
blastn <-
  #### blastn ####
  Curry(FUN = blastr:::blast, program = "blastn", task = "blastn")

#' @usage blastn_short(query, db="nt", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
#' @export
#' @rdname blastn
blastn_short <- 
  Curry(FUN = blastr:::blast, program = "blastn", task = "blastn-short")

#' @usage megablast(query, db="nt", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
#' @export
#' @rdname blastn
megablast <- 
  Curry(FUN = blastr:::blast, program = "blastn", task = "megablast")

#' @usage dc_megablast(query, db="nt", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
#' @export
#' @rdname blastn
dc_megablast <- 
  Curry(FUN = blastr:::blast, program = "blastn", task = "dc-megablast")

#' Wrapper for the NCBI Protein-Protein BLAST
#' 
#' \code{blastp} is the traditional BLASTP.
#' \code{blastp_short} is BLASTP optimised for residues shorter than 30.
#' 
#' Run \code{blastp()} without arguments to print usage and
#' arguments description.
#' 
#' @usage blastp(query, db="nr", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
#'  show_cmd=FALSE)
#'
#' @param query The sequence to search with provided as a file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param db The database to BLAST against.
#' @param out Output file for alignment. If \code{NULL} the BLAST result
#' is returned as an R character vector.
#' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
#' Other options include 6 or 7 for tabular output and 0 for the
#' traditional pairwise alignment view.
#' @param max_hits How many hits to return.
#' @param evalue Expectation value cutoff.
#' @param matrix Scoring matrix name (Default: BLOSUM62).
#' @param html Produce HTML output.
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @family blast applications
#' @export blastp blastp_short
#' @aliases blastp blastp_short
#' @examples
#' ##
blastp <- 
  #### blastp ####
  Curry(FUN = blastr:::blast, program = "blastp", task = "blastp")

#' @usage blastp_short(query, db="nr", out=NULL, outfmt=5, max_hits=20,
#'  evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
#'  show_cmd=FALSE)
#' @export
#' @rdname blastp
#' @inheritParams blastp
blastp_short <- 
  Curry(FUN = blastr:::blast, program = "blastp", task = "blastp-short")

#' Wrapper for the NCBI Translated Query-Protein Subject BLAST
#' 
#' Run \code{blastx()} without arguments to print usage and
#' arguments description.
#' 
#' @usage blastx(query, query_gencode=1, db="nr", out=NULL, outfmt=5,
#'  max_hits=20, evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, 
#'  ..., show_cmd=FALSE)
#' 
#' @param query The sequence to search with provided as a file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param query_gencode Genetic code used to translate the query
#' (default: 1).
#' @param db The database to BLAST against.
#' @param out Output file for alignment. If \code{NULL} the BLAST result
#' is returned as an R character vector.
#' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
#' Other options include 6 or 7 for tabular output and 0 for the
#' traditional pairwise alignment view.
#' @param max_hits How many hits to return.
#' @param evalue Expectation value cutoff.
#' @param matrix Scoring matrix name (Default: BLOSUM62).
#' @param html Produce HTML output.
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @family blast applications
#' @export blastx
#' @aliases blastx
#' @examples
#' ##
blastx <- 
  #### blastx ####
  Curry(FUN = blastr:::blast, program = "blastx")

#' Wrapper for the NCBI Translated Query-Protein Subject BLAST
#' 
#' Run \code{tblastx()} without arguments to print usage and
#' arguments description.
#' 
#' @usage tblastx(query, query_gencode=1, db="nr", out=NULL, outfmt=5,
#'  max_hits=20, evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE,
#'  ..., show_cmd=FALSE)
#' 
#' @param query The sequence to search with provided as a file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param query_gencode Genetic code used to translate the query
#' (default: 1).
#' @param db The database to BLAST against.
#' @param out Output file for alignment. If \code{NULL} the BLAST result
#' is returned as an R character vector.
#' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
#' Other options include 6 or 7 for tabular output and 0 for the
#' traditional pairwise alignment view.
#' @param max_hits How many hits to return.
#' @param evalue Expectation value cutoff.
#' @param matrix Scoring matrix name (Default: BLOSUM62).
#' @param html Produce HTML output.
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @family blast applications
#' @export tblastx
#' @aliases tblastx
#' @examples
#' ##
tblastx <- 
  #### tblastx ####
  Curry(FUN = blastr:::blast, program = "tblastx")

#' Wrapper for the NCBI Protein Query-Translated Subject BLAST
#' 
#' Run \code{tblastn()} without arguments to print usage and
#' arguments description.
#' 
#' @usage tblastn(query, db="nr", out=NULL, outfmt=5, max_hits=20,
#' evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
#' show_cmd=FALSE)
#' 
#' @param query The sequence to search with provided as a file name, an
#' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
#' object, or as a character vector.
#' @param db The database to BLAST against.
#' @param out Output file for alignment. If \code{NULL} the BLAST result
#' is returned as an R character vector.
#' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
#' Other options include 6 or 7 for tabular output and 0 for the
#' traditional pairwise alignment view.
#' @param max_hits How many hits to return.
#' @param evalue Expectation value cutoff.
#' @param matrix Scoring matrix name (Default: BLOSUM62).
#' @param html Produce HTML output.
#' @param remote Execute search remotely.
#' @param ... Additional parameters passed on to the BLAST commmand line
#' tools. See
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
#' for a description of common options.
#' @param show_cmd If \code{TRUE} print the constructed command line
#' instead of passing it to \code{\link{system}}.
#' 
#' @family blast applications
#' @export tblastn
#' @aliases tblastn
#' @examples
#' ##
tblastn <- 
  #### tblastn ####
  Curry(FUN = blastr:::blast, program = "tblastn")

#' Do a BLAST search using the QBLAST URL API
#' 
#' @param query An accession number, a gi, or a fasta sequence provided
#' either as an \code{\link[Biostrings]{XString}} or
#' \code{\link[Biostrings]{XStringSet}} object, or as a named character
#' vector.
#' @param program The blast application to use. One of 'megablast'
#' (default), 'blastn', 'blastp', 'rpsblast', blastx', 'tblastn', 'tblastx'.
#' @param db Blast database name. Defaults to 'nr'.
#' @param outfmt Output format. One of 'xml', 'tabular', or 'html'.
#' See Details.
#' @param max_hits The maximum number of hits to return.
#' @param evalue Expectation value threshold for saving hits.
#' @param entrez_query An Entrez query used to limit the results. 
#' @param ... The name-value pairs of parameters passed on to the QBLAST
#' URL API.
#' @param .params A named list of parameters passed on to the QBLAST
#' URL API.
#' @param display Display the query result in a web browser.
#' @param update_time How often to poll the blast server for results.
#' 
#' @return A \code{\link{blastReport-class}} object.
#'
#' @family blast applications
#' @export
#' @examples
#' ##
qblast <- function(query, program = 'megablast', db = 'nr', outfmt = 'xml',
                   max_hits = 20, entrez_query = '', evalue = 10, ...,
                   .params = list(), display = FALSE, update_time = 4) {
  if (missing(query))
    stop("No query provided")
  
  qbd <- megablast <- rpsblast <- FALSE
  program <- match.arg(program, c("megablast","blastn","blastp", "rpsblast",
                                  "blastx","tblastn", "tblastx"))
  outfmt <- match.arg(outfmt, c("xml", "tabular", "html"))
  
  if (program == "megablast") {
    program <- "blastn"
    megablast <- TRUE
  }
  
  if (program == "rpsblast") {
    program <- "blastp"
    rpsblast <- TRUE
  }
  
  if (outfmt == "xml") {
    format_type <- "XML"
    alignment_view <- "Pairwise"
    display <- display
  } 
  else if (outfmt == "tabular") {
    format_type <- "Text"
    alignment_view <- "Tabular"
    display <- display
  } 
  else if (outfmt == "html") {
    format_type <- "HTML"
    alignment_view <- "Pairwise"
    display <- TRUE
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
                             program=program,
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
  
  if (outfmt == "html") {
    displayHTML(response, unlink=TRUE)
    return(invisible(response))
  }
  else if (outfmt == "tabular") {
    if (display) displayHTML(response)
    response <- xvalue(htmlParse(response), "//pre")
    return(response)
  }
  else if (outfmt == "xml") {
    if (display) displayHTML(response, unlink=FALSE)
    return(response)
  }
}
