##' display a html file in a browser
##' 
##' @param html file path or html encoded character string
##' @param browser browser
##' @param unlink remove temporary file
##' 
##' @export 
displayHTML <- function (html, browser=getOption("browser"), unlink=TRUE)
{
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

##' Parse fasta definition lines
##' 
##' @param defline List or character vector of NCBI fasta deflines.
##' @param species Parse out species designations.
##' 
##' @importFrom stringr perl
##' @importFrom stringr str_split_fixed
##' @importFrom stringr str_extract_all
##' @importFrom stringr str_split
##' @export
parseDeflines <- function (defline, species=FALSE)
{
  # first split into identifier and description at the first blank space
  x <- str_split_fixed(unlist(defline), " ", 2)
  id <- as.list(x[,1])
  desc <- as.list(x[,2])
  
  if (species) {
    x <- str_split_fixed(unlist(desc), " \\[|\\]", 3)
    desc <- x[,1]
    sp <- x[,2]
  }
  
  ## parse identifier patterns
  ## first we extract the database tags which always are 2 or 3 lowercase
  ## letters followed by a pipe.
  db_pattern <- perl("([[:lower:]]{2,3})(?=\\|)")
  db_tag <- str_extract_all(unlist(id), db_pattern)
  db_tag[vapply(db_tag, isEmpty, logical(1))] <- "identifier"
  
  ## next we split the identifier along the database tags
  rmEmpty <- function (x) {
    if (is.atomic(x)) x <- list(x)
    Map(function (x) x[nzchar(x)], x=x)
  }
  
  str_split_list <- function (string, pattern) {  
    Map( function (string, pattern) {
      rmEmpty(str_split(string, paste(pattern, collapse="|")))[[1L]]
    }, string=string, pattern=pattern, USE.NAMES=FALSE) 
  }
  
  
  ids <- str_split_list(id, db_tag)
  id_list <- list()
  for (i in seq_along(ids)) {
    id_list[[i]] <- structure(str_split_list(ids[[i]], "\\|"), names=db_tag[[i]])
  }
  
  if (species) {
    return(list(id=id_list, desc=desc, species=sp))
  } else {
    return(list(id=id_list, desc=desc))
  }
}

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

##' Format paragraphs
##' 
##' Similar to \code{\link{strwrap}} but returns a single string with
##' linefeeds inserted
##' 
##' @param s a character vector or a list of character vectors
##' @param width a positive integer giving the column for inserting
##' linefeeds
##' @param indent an integer giving the indentation of the first line of
##' the paragraph; negative values of \code{indent} are allowed and reduce
##' the width for the first line by that value.
##' @param offset a non-negative integer giving the indentation of all
##' but the first line
##' @param split regular expression used for splitting. Defaults to
##' a whitespace character.
##' @param FORCE if \code{TRUE} words are force split if the available width
##' is too small.
##' @param FULL_FORCE Always split at the specified position.
##' 
##' @return a character vector
##' @keywords internal
linebreak <- function (s, width=getOption("width") - 2, indent=0, offset=0,
                       split=" ", FORCE=FALSE, FULL_FORCE=FALSE) {
  if (!is.character(s)) 
    s <- as.character(s)
  
  if (length(s) == 0L)
    return("")
  
  # set indent string to "" if a negative value is given
  # this lets us shrink the available width for the first line by that value
  indent_string <- blanks(ifelse(indent < 0, 0, indent))
  offset_string <- paste0("\n", blanks(offset))
  
  s <- mapply(function (s, width, offset, indent, indent_string, split, FORCE, FULL_FORCE) {
    # remove leading and trailing blanks
    # convert newlines, tabs, spaces to " "
    # find first position where 'split' applies
    if (!FULL_FORCE) {
      s <- gsub("[[:space:]]+", " ", gsub("^[[:blank:]]+|[[:blank:]]+$", "", s), perl=TRUE)
    }
    fws <- regexpr(split, s, perl=TRUE)
    if (offset + indent + nchar(s) > width) {
      # if not everything fits on one line
      if (FULL_FORCE ||
        (fws == -1 || fws >= (width - offset - indent)) && FORCE) {
        # if no whitespace or first word too long and force break
        # cut through the middle of a word
        pat1 <- paste0("^.{", width - offset - indent, "}(?=.+)")
        pat2 <- paste0("(?<=^.{", width - offset - indent, "}).+")
        leading_string <- regmatches(s, regexpr(pat1, s, perl=TRUE))
        trailing_string <- regmatches(s, regexpr(pat2, s, perl=TRUE)) 
        s <- paste0(indent_string, leading_string, offset_string,
                    linebreak(s=trailing_string, width=width, indent=0,
                              offset=offset, split=split, FORCE=FORCE,
                              FULL_FORCE=FULL_FORCE))
      } 
      else if ((fws == -1 || fws >= (width - offset + indent)) && !FORCE) {
        # if no whitespace or first word too long and NO force break
        # stop right here
        stop("Can't break in the middle of a word. Use the force!")
      }
      else {
        # break the line
        s_split <- unlist(strsplit(s, split))
        s_cum <- cumsum(nchar(s_split) + 1)
        leading_string <- 
          paste0(s_split[s_cum < width - offset - indent],
                 ifelse(split == " ", "", split), collapse=split)
        trailing_string <- 
          paste0(s_split[s_cum >= width - offset - indent], collapse=split)
        s <- paste0(indent_string, leading_string, offset_string,
                    linebreak(s=trailing_string, width=width, indent=0,
                              offset=offset, split=split, FORCE=FORCE, FULL_FORCE=FULL_FORCE))
      }
    }
    else
      # if everything fits on one line go with the string
      s
  }, s, width, offset, abs(indent), indent_string, split, FORCE, FULL_FORCE,
              SIMPLIFY=FALSE, USE.NAMES=FALSE)
  unlist(s)
}


##' create blank strings with a given number of characters
##' @seealso Examples for \code{\link{regmatches}}
##' @keywords internal
blanks <- function(n) {
  vapply(Map(rep.int, rep.int(" ", length(n)), n, USE.NAMES=FALSE),
         paste, "", collapse="")
}


##' @importFrom stringr str_pad
##' @keywords internal
wrapAln <- function (seq1,
                     ..., 
                     prefix=c(""),
                     suffix=c(""),
                     start=c(1),
                     reverse=c(FALSE),
                     sep=2)
{
  # seqs <- c(seq1, list(seq2, seq3))
  seqs <- c(seq1, list(...))
  lseqs <- vapply(seqs, nchar, numeric(1))
  
  if (!length(unique(lseqs)) == 1L)
    stop("Sequences are of different length")
  
  pref_width <- max(vapply(prefix, nchar, numeric(1))) 
  aln_start_width <- aln_end_width <-  max(c(nchar(start), nchar(unique(lseqs))))
  suf_width <- max(vapply(suffix, nchar, numeric(1)))
  offset <- pref_width + sep + aln_start_width + 1 + 1 + aln_end_width + sep + suf_width  
  
  # break up sequences  
  s <- linebreak(seqs, getOption("width") - offset, FULL_FORCE=TRUE)
  s <- str_split(s, "\n")  
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
  
  pasteAln <- function(prefix, seq_starts, s, seq_ends, suffix) {
    seq_starts[is.empty(seq_starts)] <- ""
    seq_ends[is.empty(seq_ends)] <- ""
    paste0(str_pad(prefix, pref_width, side="right"),
           blanks(sep),
           str_pad(as.character(seq_starts), aln_start_width, side="left"),
           blanks(1),
           str_pad(s, max_seq_width, side="right"),
           blanks(1),
           str_pad(as.character(seq_ends), aln_start_width, side="left"),
           blanks(sep),
           str_pad(suffix, suf_width, side="right"))
  }
  
  s <- mapply(pasteAln, prefix=prefix, seq_starts=seq_starts, s=s,
              seq_ends=seq_ends, suffix=suffix,
              SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
  s <- paste0(do.call(function (...) paste(..., sep="\n"), s),
              collapse="\n\n")
  s
}



