##' @include blast-utils.r
##' @include blast-classes.r
NULL

SysCall <- function (exec,
                     ...,
                     args = list(),
                     stdin = NULL,
                     stdout = NULL,
                     redirection = TRUE,
                     style = c("unix", "gnu"),
                     show_cmd = FALSE,
                     intern = FALSE,
                     input = NULL)
{  
  isFALSE <- function (x) identical(FALSE, x)
  
  args <- merge(list(...), args)
  style <- match.arg(style)
  
  if (is.null(stdin)) {
    stdin <- ""
  }
  else if (!is.null(stdin) && redirection) {
    stdin <- paste("<", stdin)
  }
  
  if (is.null(stdout)) {
    stdout <- ""
  }
  else {
    stdout <- paste(">", stdout)
  }
   
  args[vapply(args, isTRUE, logical(1))] <- ""
  args[vapply(args, isFALSE, logical(1))] <- NULL
  args[vapply(args, is.null, logical(1))] <- NULL
  args <- switch(style,
                 unix=paste(str_trim(sprintf("-%s %s", names(args), args)),
                            collapse=" "),
                 gnu=paste(str_trim(sprintf("--%s %s", names(args), args)),
                           collapse=" ")
                 )
  if (show_cmd)
    print(str_trim(paste(exec, args, stdin, stdout)))
  else
    return(system(str_trim(paste(exec, args, stdin, stdout)),
                  intern = intern, input = input))
}

Curry <- function (FUN, ...)
{
  args <- match.call(expand.dots=FALSE)$...
  args$... <- as.name("...")
  
  env <- new.env(parent=parent.frame())
  
  if (is.name(FUN)) {
    fname <- FUN
  } else if (is.character(FUN)) {
    fname <- as.name(FUN)
  } else if (is.function(FUN)) {
    fname <- as.name("FUN")
    env$FUN <- FUN
  } else {
    stop("FUN not function or name of function")
  }
  
  curry_call <- as.call(c(list(fname), args))
  
  f <- eval(call("function", as.pairlist(alist(... = )), curry_call))
  environment(f) <- env
  f
}

##' Wrapper for NCBI makeblastdb
##' 
##' @param input_file Input file/database name.
##' @param input_type Type of data specified in input file.
##' @param dbtype Molecule type of target db.
##' @param ... further arguments passed to makeblastdb.
##' @param show_log print log file.
##' @param show_cmd print the command line instead of executing it.
##' 
##' @family blast applications
##' 
##' @export
makeblasttdb <- function(input_file,
                         input_type=c("fasta","blastdb","asn1_bin","asn1_txt"),
                         dbtype=c("nucl","prot"),
                         ...,
                         show_log=TRUE,
                         show_cmd=FALSE)
{
  if (missing(input_file))
    return(SysCall("makeblastdb", help=TRUE, redirection=FALSE))
  
  if (!file.exists(input_file))
    stop(sprintf("File %s does not exist", sQuote(input_file)))

  input_type <- match.arg(input_type)
  dbtype <- match.arg(dbtype)
  
  o <- list(...)
  if (!is.null(o$logfile)) logfile <- o$logfile else logfile <- ""
  
  logfile <- "logfile"
  
  SysCall(exec="makeblastdb", infile=NULL, outfile=NULL,
          `in`=input_file, input_type=input_type,
          dbtype=dbtype, ..., style="unix", show_cmd=show_cmd)
  
  if (nzchar(logfile) && isTRUE(show_log))
    cat(paste(readLines(logfile), collapse="\n"))
}

##' Wrapper for the new NCBI BLAST+ tools
##' 
##' @param program One of blastn, blastp, blastx, tblastn, or tblastx
##' @param query Query provided as an input file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param db Blast database name (defaults to 'nt' for blastn and 'nr' for
##' the rest).
##' @param outfmt Output format (defaults to XML output).
##' @param max_hits Maximum number of hits  to return. Defaults to 20.
##' Sets the '-max_target_seqs' parameter internally.
##' @param strand Query strand(s) to seach against database.
##' @param ... Arguments passed on to the blast commmand line tools.
##' @param intern Set \code{TRUE} if no '-out' argument is specified.
##' Captures the blast output in an R character vector.
##' @param input Used to pass a character vector to the standard input.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @keywords internal
blast <- function (program = c("blastn", "blastp", "blastx", "tblastn", "tblastx"),
                   query,
                   db,
                   outfmt,
                   max_hits,
                   strand = c("both", "plus", "minus"),
                   ...,
                   intern = FALSE,
                   input = NULL,
                   show_cmd = FALSE)
{  
  program <- match.arg(program)
  strand <- match.arg(strand)
  
  ### dealing with query
  if (missing(query))
    return(SysCall(program, help=TRUE, intern=FALSE))
  
  if (is(query, "XStringSet")) {
    input <- paste0(">", names(query[1L]), "\n", toString(query[[1L]]))
    query <- NULL
  }
  else if (is(query, "XString")) {
    input <- paste0(">", names(query[1L]), "\n", toString(query))
    query <- NULL
  }
  else if (is.vector(query)) {
    if (!file.exists(query)) {
      input <- as.character(query)
      query <- NULL
    }
    else {
      input <- NULL
      query <- query
    }
  } 
  else {
    stop(sprintf("Objects of class %s are not supported for query",
                 sQuote(class(query))))
  }
  
  ## set a number of defaults different from the internal defaults of
  ## the blast applications
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
  }
  else if (as.integer(outfmt) >= 5) {
    num_descriptions <- NULL
    num_alignments <- NULL
    max_target_seqs <- max_hits
  }
  
  args <- merge(list(...),
               # args <-
                list(query=query,
                     db=db,
                     outfmt=outfmt,
                     num_descriptions=num_descriptions,
                     num_alignments=num_alignments,
                     max_target_seqs=max_target_seqs,
                     strand=strand))
  
  # remove stdin', 'stdout' from the arguments list
  stdin <- args[["stdin"]]
  args[["stdin"]] <- NULL
  stdout <- args[["stdout"]]
  args[["stdout"]] <- NULL
  
  # check if 'out' is specified, otherwise return results internally
  intern <- if (is.null(args[["out"]])) TRUE else FALSE
  
  res <- SysCall(exec=program, args=args,
                 stdin=stdin, stdout=stdout,
                 redirection=if (is.null(stdin) && is.null(stdout)) FALSE else TRUE,
                 style="unix", show_cmd=show_cmd,
                 intern=intern, input=input)
  
  res <- paste(res, collapse="\n")
  return(res)
}

##' Wrapper for the NCBI Nucleotide-Nucleotide BLAST
##' 
##' \itemize{
##' \item{\code{blastn}} is the traditional BLASTN requiring an exact
##' match of 11.
##' \item{\code{blastn_short}} is BLASTN optimised for sequences shorter
##' than 50 bases.
##' \item{\code{megablast}} is the traditional megablast used to find
##' very similar sequences.
##' \item{\code{dc_megablast}} is discontiguous megablast used to find
##' more distant (e.g. interspecies) sequences.
##' }
##' 
##' Run \code{blastn()} without arguments to print usage and
##' arguments description.
##' 
##' @usage blastn(query, db="nt", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
##'
##' @param query The sequence to search with provided as a file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param db The database to BLAST against.
##' @param out Output file for alignment. If \code{NULL} the BLAST result
##' is returned as an R character vector.
##' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
##' Other options include 6 or 7 for tabular output and 0 for the
##' traditional pairwise alignment view.
##' @param max_hits How many hits to return.
##' @param evalue Expectation value cutoff.
##' @param html Produce HTML output.
##' @param remote Execute search remotely.
##' @param ... Additional parameters passed on to the BLAST commmand line
##' tools. See
##' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
##' for a description of common options.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @family blast applications
##' @export blastn blastn_short megablast dc_megablast
##' @aliases blastn blastn_short megablast dc_megablast
##' @examples
##' ##
blastn <-
  #### blastn ####
  Curry(FUN = blast, program = "blastn", task = "blastn")

##' @usage blastn_short(query, db="nt", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
##' @export
##' @rdname blastn
blastn_short <- Curry(FUN = blast, program = "blastn", task = "blastn-short")

##' @usage megablast(query, db="nt", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
##' @export
##' @rdname blastn
megablast <- Curry(FUN = blast, program = "blastn", task = "megablast")

##' @usage dc_megablast(query, db="nt", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, html=FALSE, remote=FALSE, ..., show_cmd=FALSE)
##' @export
##' @rdname blastn
dc_megablast <- Curry(FUN = blast, program = "blastn", task = "dc-megablast")

##' Wrapper for the NCBI Protein-Protein BLAST
##' 
##' \code{blastp} is the traditional BLASTP.
##' \code{blastp_short} is BLASTP optimised for residues shorter than 30.
##' 
##' Run \code{blastp()} without arguments to print usage and
##' arguments description.
##' 
##' @usage blastp(query, db="nr", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
##'  show_cmd=FALSE)
##'
##' @param query The sequence to search with provided as a file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param db The database to BLAST against.
##' @param out Output file for alignment. If \code{NULL} the BLAST result
##' is returned as an R character vector.
##' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
##' Other options include 6 or 7 for tabular output and 0 for the
##' traditional pairwise alignment view.
##' @param max_hits How many hits to return.
##' @param evalue Expectation value cutoff.
##' @param matrix Scoring matrix name (Default: BLOSUM62).
##' @param html Produce HTML output.
##' @param remote Execute search remotely.
##' @param ... Additional parameters passed on to the BLAST commmand line
##' tools.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @family blast applications
##' @export blastp blastp_short
##' @aliases blastp blastp_short
##' @examples
##' ##
blastp <- 
  #### blastp ####
  Curry(FUN = blast, program = "blastp", task = "blastp")

##' @usage blastp_short(query, db="nr", out=NULL, outfmt=5, max_hits=20,
##'  evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
##'  show_cmd=FALSE)
##' @export
##' @rdname blastp
##' @inheritParams blastp
blastp_short <- Curry(FUN = blast, program = "blastp", task = "blastp-short")

##' Wrapper for the NCBI Translated Query-Protein Subject BLAST
##' 
##' Run \code{blastx()} without arguments to print usage and
##' arguments description.
##' 
##' @usage blastx(query, query_gencode=1, db="nr", out=NULL, outfmt=5,
##'  max_hits=20, evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, 
##'  ..., show_cmd=FALSE)
##' 
##' @param query The sequence to search with provided as a file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param query_gencode Genetic code used to translate the query
##' (default: 1).
##' @param db The database to BLAST against.
##' @param out Output file for alignment. If \code{NULL} the BLAST result
##' is returned as an R character vector.
##' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
##' Other options include 6 or 7 for tabular output and 0 for the
##' traditional pairwise alignment view.
##' @param max_hits How many hits to return.
##' @param evalue Expectation value cutoff.
##' @param matrix Scoring matrix name (Default: BLOSUM62).
##' @param html Produce HTML output.
##' @param remote Execute search remotely.
##' @param ... Additional parameters passed on to the BLAST commmand line
##' tools. See
##' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
##' for a description of common options.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @family blast applications
##' @export blastx
##' @aliases blastx
##' @examples
##' ##
blastx <- 
  #### blastx ####
  Curry(FUN = blast, program = "blastx")

##' Wrapper for the NCBI Translated Query-Protein Subject BLAST
##' 
##' Run \code{tblastx()} without arguments to print usage and
##' arguments description.
##' 
##' @usage tblastx(query, query_gencode=1, db="nr", out=NULL, outfmt=5,
##'  max_hits=20, evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE,
##'  ..., show_cmd=FALSE)
##' 
##' @param query The sequence to search with provided as a file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param query_gencode Genetic code used to translate the query
##' (default: 1).
##' @param db The database to BLAST against.
##' @param out Output file for alignment. If \code{NULL} the BLAST result
##' is returned as an R character vector.
##' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
##' Other options include 6 or 7 for tabular output and 0 for the
##' traditional pairwise alignment view.
##' @param max_hits How many hits to return.
##' @param evalue Expectation value cutoff.
##' @param matrix Scoring matrix name (Default: BLOSUM62).
##' @param html Produce HTML output.
##' @param remote Execute search remotely.
##' @param ... Additional parameters passed on to the BLAST commmand line
##' tools. See
##' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
##' for a description of common options.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @family blast applications
##' @export tblastx
##' @aliases tblastx
##' @examples
##' ##
tblastx <- 
  #### tblastx ####
  Curry(FUN = blast, program = "tblastx")

##' Wrapper for the NCBI Protein Query-Translated Subject BLAST
##' 
##' Run \code{tblastn()} without arguments to print usage and
##' arguments description.
##' 
##' @usage tblastn(query, db="nr", out=NULL, outfmt=5, max_hits=20,
##' evalue=10, matrix="BLOSUM62", html=FALSE, remote=FALSE, ...,
##' show_cmd=FALSE)
##' 
##' @param query The sequence to search with provided as a file name, an
##' \code{\link[Biostrings]{XString}} or \code{\link[Biostrings]{XStringSet}}
##' object, or as a character vector.
##' @param db The database to BLAST against.
##' @param out Output file for alignment. If \code{NULL} the BLAST result
##' is returned as an R character vector.
##' @param outfmt Output format, Integer 1-11. Default is 5 for XML output.
##' Other options include 6 or 7 for tabular output and 0 for the
##' traditional pairwise alignment view.
##' @param max_hits How many hits to return.
##' @param evalue Expectation value cutoff.
##' @param matrix Scoring matrix name (Default: BLOSUM62).
##' @param html Produce HTML output.
##' @param remote Execute search remotely.
##' @param ... Additional parameters passed on to the BLAST commmand line
##' tools. See
##' \href{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.4_User_manual}{here}
##' for a description of common options.
##' @param show_cmd If \code{TRUE} print the constructed command line
##' instead of passing it to \code{\link{system}}.
##' 
##' @family blast applications
##' @export tblastn
##' @aliases tblastn
##' @examples
##' ##
tblastn <- 
  #### tblastn ####
  Curry(FUN = blast, program = "tblastn")

##' Do a BLAST search using the QBLAST URL API
##' 
##' @param query An accession number, a gi, or a fasta sequence provided
##' either as an \code{\link[Biostrings]{XString}} or
##' \code{\link[Biostrings]{XStringSet}} object, or as a named character
##' vector.
##' @param program The blast application to use. One of 'megablast'
##' (default), 'blastn', 'blastp', 'rpsblast', blastx', 'tblastn', 'tblastx'.
##' @param db Blast database name. Defaults to 'nr'.
##' @param outfmt Output format. One of 'xml', 'tabular', or 'html'.
##' See Details.
##' @param max_hits The maximum number of hits to return.
##' @param evalue Expectation value threshold for saving hits.
##' @param entrez_query An Entrez query used to limit the results. 
##' @param ... The name-value pairs of parameters passed on to the QBLAST
##' URL API.
##' @param .params A named list of parameters passed on to the QBLAST
##' URL API.
##' @param display Display the query result in a web browser.
##' @param update_time How often to poll the blast server for results.
##' 
##' @return A \code{\link{blastReport-class}} object.
##'
##' @importFrom stringr str_match
##'
##' @family blast applications
##' @export
##' @examples
##' ##
qblast <- function(query,
                   program=c("megablast","blastn","blastp",
                             "rpsblast","blastx","tblastn",
                             "tblastx"),
                   db="nr",
                   outfmt=c("xml", "tabular", "html"),
                   max_hits=20,
                   entrez_query="",
                   evalue=10,
                   ...,
                   .params=list(),
                   display=FALSE,
                   update_time=4)
{
  if (missing(query))
    stop("No query provided")
  
  qbd <- megablast <- rpsblast <- FALSE
  program <- match.arg(program)
  outfmt <- match.arg(outfmt)
  
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
  } 
  else if (!is.null(names(query))) {
    query <- paste0(">", names(query)[[1L]], "\n", as.character(query[[1L]]))
    qbd <- TRUE
  } 
  else {
    query <- as.character(query)
  }
  
  if (nchar(entrez_query) > 0L) {
    .params <- merge(.params, list(entrez_query=.escape(entrez_query)))
  }
  
  .params <- merge(list(query=query,
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
  xpath <- "/descendant::comment()[contains(.,'RID =')]"
  r <- xpathSApply(htmlParse(response), xpath, xmlValue)
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
    
    xpath <- "/descendant::comment()[contains(.,'Status=')]"
    status <- xpathSApply(htmlParse(response), xpath, xmlValue)
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
      xpath <- "/descendant::comment()[contains(.,'ThereAreHits=')]"
      there_are_hits <- xpathSApply(htmlParse(response), xpath, xmlValue)
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
    response <- xpathApply(htmlParse(response), "//pre", xmlValue)[[1L]]
    return(response)
  }
  else if (outfmt == "xml") {
    if (display) displayHTML(response, unlink=FALSE)
    return(response)
  }
}
