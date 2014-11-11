#' @include utils.r
#' @importFrom stringr str_split_fixed str_extract_all str_detect str_match perl
NULL

## NCBI BLAST defline database tags
.tags  <- c("lcl","gb","emb","pir","sp","ref","gnl", "gi","dbj","prf","pdb","pat","bbs")

#' @keywords internal
.escape <- function(s, httpPOST=FALSE) {
  if (httpPOST) {
    s <- gsub("\\s+", " ", s)
    s <- gsub("+", " ", s, fixed=TRUE)
  } else {
    s <- gsub("\\s+", "\\+", s)
  }
  s <- paste(strsplit(s, '\"', fixed=TRUE)[[1L]], collapse="%22")
  s <- gsub(">", "%3E", s)
  s <- gsub("\\n", "%0D%0A", s)
  s <- gsub("\\|", "%7C", s)
  s <- gsub("\\#", "%23", s)
  s <- gsub("\\+(and)\\+|\\+(or)\\+|\\+(not)\\+","\\+\\U\\1\\U\\2\\U\\3\\+", s, perl=TRUE)
  s
}

#' Construct deflines
#' @keywords internal
make_deflines <- function(query, prefix = "lcl") {
  ##
  if (class(query) %in% c("gbFeatureList", "gbFeature")) {
    check_biofiles() ## conditionally load biofiles
    id <- paste0(prefix, "|", index(query))
    desc <- paste0(unlist(qualif(query, "locus_tag")), " [", unlist(qualif(query, "product")), "]")
    parse_defline <- TRUE
  } else if (is(query, "XStringSet")) {
    # test if the XStrings follow the naming convention
    # from biofiles: accn.key.idx
    p <- "[[:alnum:]]+\\.[[:alnum:]]+\\.[[:digit:]]+"
    n <- names(query)
    if (!is.null(n) && all(grepl(p, n))) {
      sp <- vapply(n, function (x) strsplit(x, "\\.")[[1L]], character(3))
      id <- paste0(prefix, "|", sp[3L, ])
      desc <- paste0(sp[1L, ], " [", sp[2L, ], "]")
      parse_defline <- TRUE
    } else {
      id <- names(query)
      desc <- NULL
      parse_defline <- FALSE
    }
  } else {
    id <- names(query)
    desc <- NULL
    parse_defline <- FALSE
  }
  if (is.null(id)) {
    id <- paste0("Query_", seq_along(query))
  }
  list(defline=paste(id, desc), parse_defline=parse_defline)
}

#' @importFrom Biostrings fasta.info
#' @keywords internal
make_blast_query <- function(query, transl = FALSE) {
  ## Set up tempfile to use as input for blast
  tmp <- tempfile(fileext=".fa")
  if (is.string(query) && tryCatch(is.readable(query), assertError = function(e) FALSE )) {
    # query must be the path to a valid FASTA file
    nqueries <- length(fasta.info(query))
    # copy query to tempfile because blast deletes the query file
    # after it's finished
    file.copy(query, tmp)
    return(list(query=tmp, nqueries=nqueries, parse_defline=FALSE))
  }
  if (class(query) %in% c("gbReportList","gbReport","gbFeatureList","gbFeature")) {
    check_biofiles() ## conditionally load biofiles
    seq <- getSequence(query)
  } else if (inherits(query, "XStringSet") || inherits(query, "XString")) {
    seq <- as(query, "XStringSet")
  } else if (is.vector(query) && is.character(query)) {
    seq <- query
  } else {
    stop("Objects of class ", sQuote(class(query)), " are not supported as query.")
  }
  seqnames <- make_deflines(seq)
  writeLines(paste0(paste0(">", seqnames$defline, "\n", as.character(seq)), collapse="\n"), tmp)
  list(query=tmp, nqueries=length(seqnames$defline), parse_deflines=seqnames$parse_defline) 
}

#' @keywords internal
wrap_alignment <- function(seq1, ...,  prefix = "", suffix = "",
                          start = 1, reverse = FALSE, sep = 2) {
  
  paste_alignment <- function(prefix, seq_starts, s, seq_ends, suffix) {
    paste0(
      pad(prefix, pref_width, "right"), blanks(sep),
      pad(seq_starts, aln_start_width, "left"), blanks(1),
      pad(s, max_seq_width, "right"), blanks(1),
      pad(seq_ends, aln_start_width, "left"), blanks(sep),
      pad(suffix, suf_width, "right")
    )
  }
  
  # seqs <- c(seq1, list(seq2, seq3))
  seqs <- c(list(seq1), list(...))
  lseqs <- vapply(seqs, nchar, numeric(1))
  
  if (!length(unique(lseqs)) == 1L) {
    stop("Sequences are of different length")
  }
  
  pref_width <- max(vapply(prefix, nchar, numeric(1))) 
  aln_start_width <- aln_end_width <-  max(c(nchar(start), nchar(unique(lseqs))))
  suf_width <- max(vapply(suffix, nchar, numeric(1)))
  offset <- pref_width + sep + aln_start_width + 1 + 1 + aln_end_width + sep + suf_width  
  
  # break up sequences  
  s <- linebreak(seqs, getOption("width") - offset - 2, FULL_FORCE = TRUE)
  s <- strsplit(s, "\n")  
  seq_widths <- nchar(s[[1L]])
  max_seq_width <- max(seq_widths)
  
  seq_starts <- mapply(function(start, rev) {
    x <- Reduce("+", seq_widths, init = start, right = rev, accumulate = TRUE)
    x <- x[-which.max(x)]
    x
  }, start = start, rev = reverse, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  new_starts <- mapply(function(s, rev) if (rev) s[length(s) - 1] - 1 else s[2] - 1,
                       s = seq_starts, rev = reverse)
  
  seq_ends <- mapply(function(start, rev) {
    x <- Reduce("+", seq_widths, init = start, right = rev, accumulate = TRUE)
    x <- x[-which.max(x)]
  }, start = new_starts, rev = reverse, SIMPLIFY = FALSE, USE.NAMES = TRUE)  
  
  tmp <- seq_ends[reverse]
  seq_ends[reverse] <- seq_starts[reverse]
  seq_starts[reverse] <- tmp
 
  seq_starts[vapply(seq_starts, function(x) length(x) == 0, FALSE)] <- ""
  seq_ends[vapply(seq_ends, function(x) length(x) == 0, FALSE)] <- ""

  s <- .mapply(paste_alignment, list(prefix = prefix, seq_starts = seq_starts,
                                     s = s, seq_ends = seq_ends, suffix = suffix),
               NULL)
  paste0(do.call(function(...) paste(..., sep="\n"), s), collapse="\n\n")
}

