
# blast-parser -----------------------------------------------------------

#' @include blast-utils.r
#' @include blast-classes.r
NULL

#' Parse xml blast output
#' 
#' @param blast_output Blast output in XML format
#' (File or character vector) 
#' 
#' @return A \code{\link{blastReport-class}} object.
#' 
#' @importFrom Biostrings BStringSet
#' @importFrom Biostrings BString
#' @export
parseBlastXml <- function (blast_output) {
  res <- list()
#   blout <- blast_output[[1]]
  for (blout in blast_output) {
    doc <- xmlRoot(xmlParse(blout))
    program <- xvalue(doc, '//BlastOutput_program')
    version <- xvalue(doc, '//BlastOutput_version')
    reference <- xvalue(doc, '//BlastOutput_reference')
    db <- xvalue(doc, '//BlastOutput_db')
    iter_num <- xvalue(doc, '//BlastOutput_iter-num', as="integer")
    message <- xvalue(doc, '//BlastOutput_message')
    
    # BlastOutput_query
    qid <- parseDeflines(xvalue(doc, '//Iteration_query-ID'))[["id"]]
    qdef <- xvalue(doc, '//Iteration_query-def')
    qlen <- xvalue(doc, '//Iteration_query-len', as='integer')
    qseq <- xvalue(doc, '//Iteration_query-seq') # optional
    query <- merge_list(qid[[1L]], list(
      def=qdef,
      len=qlen,
      seq=if (not.na(qseq)) BString(qseq) else NULL
    ))
    
    # BlastOutput/BlastOutput_param/Parameters
    params <- as.list(setNames(xvalue(doc, '//Parameters/*'),
                               gsub("-", "_", strsplitN(xname(doc, '//Parameters/*'),
                                                        "_", 2))))
    
    ## BlastOutput/BlastOutput_Iterations//Statistics
    stats <- as.list(setNames(xvalue(doc, '//Statistics/*', as='numeric'),
                              gsub("-", "_", strsplitN(xname(doc, '//Statistics/*'),
                                                       "_", 2))))
    
    ## Hits
    hits <- getNodeSet(doc, "//Hit")
    hit_list <- list()
    for (i in seq_along(hits)) {
      hit <- xmlDoc(hits[[i]]) 
      ## parse HSPs
      hsp_obj <- .hsp(num = xvalue(hit, '//Hsp_num', as='integer'),
                      bit_score = xvalue(hit, '//Hsp_bit-score', as='numeric'),
                      score = xvalue(hit, "//Hsp_score", as='integer'),
                      evalue = xvalue(hit, "//Hsp_evalue",, as='numeric'),
                      query_from = xvalue(hit, "//Hsp_query-from", as='integer'),
                      query_to = xvalue(hit, "//Hsp_query-to", as='integer'),
                      hit_from = xvalue(hit, "//Hsp_hit-from", as='integer'),
                      hit_to = xvalue(hit, "//Hsp_hit-to", as='integer'),
                      pattern_from = xvalue(hit, "//Hsp_pattern-from", as='integer'),
                      pattern_to = xvalue(hit, "//Hsp_pattern-to", as='integer'),
                      query_frame = xvalue(hit, "//Hsp_query-frame", as='integer'),
                      hit_frame = xvalue(hit, "//Hsp_hit-frame", as='integer'),
                      identity = xvalue(hit, "//Hsp_identity", as='integer'),
                      positive = xvalue(hit, "//Hsp_positive", as='integer'),
                      gaps = xvalue(hit, "//Hsp_gaps", as='integer'),
                      density = xvalue(hit, "//Hsp_density", as='numeric'),
                      align_len = xvalue(hit, "//Hsp_align-len", as='integer'),
                      qseq = {
                        qseq <- BStringSet(xvalue(hit, "//Hsp_qseq"))
                        names(qseq) <- paste0("hsp", xvalue(hit, "//Hsp_num"))
                        qseq }, 
                      hseq = {
                        hseq <- BStringSet(xvalue(hit, "//Hsp_hseq"))
                        names(hseq) <- paste0("hsp", xvalue(hit, "//Hsp_num"))
                        hseq },
                      midline = xvalue(hit, "//Hsp_midline"),
                      percent_identity = xvalue(hit, "//Hsp_percent-identity", as='numeric'),
                      mismatch_count = xvalue(hit, "//Hsp_mismatch-count", as='integer'))
      
      ## parse hits
      id <- paste(xvalue(hit, "//Hit_id"), xvalue(hit, "//Hit_def"))
      id <- parseDeflines(defline=strsplit(id, " >")[[1L]])
      hit_obj <- .hit(num = xvalue(hit, "//Hit_num", as='integer'),
                      id = id$id,
                      desc = id$desc,
                      accn = xvalue(hit, "//Hit_accession"),
                      len = xvalue(hit, "//Hit_len",as='integer'),
                      hsp = hsp_obj)
      
      hit_list[[i]] <- hit_obj
    }
    
    
    res <- c(res, .blastReport(program=program, version=version,
                               reference=reference, db=db, query=query,
                               iter_num=iter_num, hits=hit_list, params=params,
                               stats=stats, message=message, data=blout))  
  }
  
  return(res)
}

#' read a blast hit table
#' 
#' @param blast_output Blast output in XML format
#' (File or character vector)
#' 
#' @return A \code{\link{blastReport-class}} object.
#' 
#' @export
parseBlastTabular <- function (blast_output) {
  res <- list()
  for (blout in blast_output) {
    if (is.character(blout) && file.exists(blout)) {
    cf <- count.fields(blout, sep="\t", comment.char="#")
    file_path <- file(blout, open="r")
  } else {
    cf <- count.fields(textConnection(as.character(blout)),
                       sep="\t", comment.char="#")
    file_path <- textConnection(as.character(blout))
  }

  if (!all(cf == 12)) {
    stop(sprintf("Malformed blast table. %s columns in row %s.\n",
                 sQuote(cf[cf > 12]), sQuote(which(cf > 12))))
  }
  
  hasComment <- TRUE
  comStr <- NULL
  while (hasComment) {
    line <- readLines(file_path, n=1)
    if (hasComment <- str_detect(line, "^#")) {
      comStr <- c(comStr, str_match(line, "[^(# *)].*"))
    }
  }
  
  pushBack(line, connection=file_path)  
  hit_table <- 
    read.table(file_path, header=FALSE, sep="\t",
               as.is=TRUE, nrows=length(cf),
               col.names=c("qid", "sid", "pident",
                           "length", "mismatch", "gapopen",
                           "qstart", "qend", "sstart",
                           "send","evalue","bitscore"),
               colClasses=c("character", "character", "numeric",
                            "integer", "integer", "integer",
                            "integer", "integer", "integer",
                            "integer", "numeric", "numeric"))

  ## parse subjectIds in 'hit_table'
  all_ids <- with(hit_table, strsplit(sid, "\\|"))
  gi  <- vapply(all_ids, '[', 2, FUN.VALUE=character(1))
  source_tag <- vapply(all_ids, '[', 3, FUN.VALUE=character(1))
  accn <- ifelse(grepl("pdb", source_tag),
                 paste0(sapply(all_ids, '[', 4), "_", sapply(all_ids, '[', 5)),
                 sapply(all_ids, '[', 4))
  if (all(is.na(accn)))
    accn <- hit_table[["sid"]]
  
  neg_log_evalue <- with(hit_table, -log(as.numeric(evalue)))
  neg_log_evalue[is.infinite(neg_log_evalue)] <- -log(1e-323)
  
  if (not.null(comStr) && length(comStr) == 5) {
    program <- comStr[1]
    query <- str_split_fixed(comStr[2], "Query: ", 2)[,2]
    db <- str_split_fixed(comStr[3], "Database: ", 2)[,2]
  } else {
    program <- query <- db <- NA_character_
  }
    
    close(file_path)
    
    res <- c(res, .blastTable(program = program, db = db, query = query,
                              bitscore = as.numeric(hit_table[["bitscore"]]),
                              evalue = as.numeric(hit_table[["evalue"]]),
                              mlog.evalue = neg_log_evalue,
                              gi = gi,
                              accession = accn,
                              table = hit_table))
  }
  
  return(res)
}
