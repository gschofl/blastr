#' @importFrom biofiles sequence
#' @importFrom biofiles strand
#' @importFrom biofiles index
#' @importFrom biofiles qualif
#' @importFrom rmisc linebreak
#' @importFrom rmisc pad
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_detect
#' @importFrom stringr str_match
#' @importFrom stringr perl
#' @importFrom Biostrings translate
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings BStringSet
NULL


#' display a html file in a browser
#' 
#' @param html file path or html encoded character string
#' @param browser browser
#' @param unlink remove temporary file
#' 
#' @export 
displayHTML <- function (html, browser=getOption("browser"), unlink=TRUE) {
  if (!file.exists(html)) { 
    f_tmp <- tempfile(fileext=".html")
    writeLines(html, f_tmp)
  } else {
    f_tmp <- html
    unlink <- FALSE
  }
  
  browseURL(url=f_tmp, browser=browser)
  
  if (unlink) {
    Sys.sleep(2)
    unlink(f_tmp)
  }
}


#' @keywords internal
.escape <- function (s, httpPOST=FALSE) {
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
#' deparse NCBI fasta definition lines
#' 
#' @param ids List of identifiers.
#' @param descs List of description lines.
#' 
#' @keywords internal
deparseDeflines <- function(ids, descs) {

  deflines <- list()
  for (i in seq_along(ids)) {
    deflines[[i]] <-
      paste0(
        paste0(names(ids[[i]]), "|", paste0(ids[[i]][[1]], collapse="|") , collapse="|"),
        "| ",
        descs[[i]])
  }
  
  return(deflines)
}


#' Construct deflines
#' @keywords internal
make_deflines <- function (x, prefix = "lcl") {
  if (is(x, "gbFeature") || is(x, "gbFeatureList")) {
    id <- paste0(prefix, "|", index(x))
    desc <- paste0(unlist(qualif(x, "locus_tag")),
                   " [", unlist(qualif(x, "product")), "]")
    parse_defline <- TRUE
  } else if (is(x, "XStringSet")) {
    # test if the XStrings follow the naming convention
    # from biofiles: accn.key.idx
    p <- "[[:alnum:]]+\\.[[:alnum:]]+\\.[[:digit:]]+"
    n <- names(x)
    if (!is.null(n) && all(grepl(p, n))) {
      sp <- vapply(n, function (x) strsplit(x, "\\.")[[1L]], character(3))
      id <- paste0(prefix, "|", sp[3L, ])
      desc <- paste0(sp[1L, ], " [", sp[2L, ], "]")
      parse_defline <- TRUE
    } else {
      id <- names(x)
      desc <- NULL
      parse_defline <- FALSE
    }
  } else {
    id <- names(x)
    desc <- NULL
    parse_defline <- FALSE
  }
  if (is.null(id)) {
    id <- paste0("Query_", seq_along(x))
  }
  list(defline=paste(id, desc), parse_defline=parse_defline)
}


#' @keywords internal
make_blast_query <- function (x, transl = FALSE) {
  ## when x is a path to a FASTA file
  if (is.string(x) && is.readable(x)) {
    return( list(query=x, input=NULL, parse_defline=FALSE) )
  }
  if (is(x, "gbFeatureList") || is(x, "gbFeature")) {
    seq <- biofiles::sequence(x)
  } else if (is(x, "XString") || is(x, "XStringSet")) {
    seq <- as(x, "XStringSet")
  } else if (is.vector(x) && is.character(x)) {
    seq <- x
  } else {
    stop(sprintf("Objects of class %s are not supported as query",
                 sQuote(class(query))))
  }
  
  seqnames <- make_deflines(x, prefix="lcl")
  
  if (transl && class(x) %in% c("gbFeature","gbFeatureList")) {
    plus <- biofiles::strand(x) == 1
    minus <- !plus
    seq <- c(translate(reverseComplement(seq[minus])), translate(seq[plus]))
    seqnames$defline <- c(seqnames$defline[minus], seqnames$defline[plus])
  }
  
  input <- paste0(paste0(">", seqnames$defline, "\n", as.character(seq)),
                  collapse="\n")
  list(query=NULL, input=input, deflines=seqnames$defline,
       parse_deflines=seqnames$parse_defline) 
}


#' @keywords internal
wrapAlignment <- function (seq1, ...,  prefix=c(""), suffix=c(""),
                           start=c(1), reverse=c(FALSE), sep=2) {
  # seqs <- c(seq1, list(seq2, seq3))
  seqs <- c(list(seq1), list(...))
  lseqs <- vapply(seqs, nchar, FUN.VALUE=numeric(1))
  
  if (!length(unique(lseqs)) == 1L)
    stop("Sequences are of different length")
  
  pref_width <- max(vapply(prefix, nchar, numeric(1))) 
  aln_start_width <- aln_end_width <-  max(c(nchar(start), nchar(unique(lseqs))))
  suf_width <- max(vapply(suffix, nchar, numeric(1)))
  offset <- pref_width + sep + aln_start_width + 1 + 1 + aln_end_width + sep + suf_width  
  
  # break up sequences  
  s <- linebreak(seqs, getOption("width") - offset - 2, FULL_FORCE=TRUE)
  s <- strsplit(s, "\n")  
  seq_widths <- nchar(s[[1L]])
  max_seq_width <- max(seq_widths)
  
  seq_starts <-
    mapply(function (start, rev) {
      x <- Reduce("+", seq_widths, init=start, right=rev, accumulate=TRUE)
      x <- x[-which.max(x)]
      x
    }, start=start, rev=reverse, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
  new_starts <- 
    mapply( function (s, rev) if (rev) s[length(s) - 1] - 1 else s[2] - 1,
            s=seq_starts, rev=reverse)
  
  seq_ends <-
    mapply( function (start, rev) {
      x <- Reduce("+", seq_widths, init=start, right=rev, accumulate=TRUE)
      x <- x[-which.max(x)]
    }, start=new_starts, rev=reverse, SIMPLIFY=FALSE, USE.NAMES=FALSE)  
  
  tmp <- seq_ends[reverse]
  seq_ends[reverse] <- seq_starts[reverse]
  seq_starts[reverse] <- tmp
  
  pasteAlignment <- function(prefix, seq_starts, s, seq_ends, suffix) {
    seq_starts[is_empty(seq_starts)] <- ""
    seq_ends[is_empty(seq_ends)] <- ""
    paste0(pad(prefix, pref_width, "right"), blanks(sep),
           pad(seq_starts, aln_start_width, "left"), blanks(1),
           pad(s, max_seq_width, "right"), blanks(1),
           pad(seq_ends, aln_start_width, "left"), blanks(sep),
           pad(suffix, suf_width, "right"))
  }
  
  s <- mapply(pasteAlignment, prefix=prefix, seq_starts=seq_starts, s=s,
              seq_ends=seq_ends, suffix=suffix,
              SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
  s <- paste0(do.call(function (...) paste(..., sep="\n"), s),
              collapse="\n\n")
  s
}

