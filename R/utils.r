#' @import methods
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.readable
#' @importFrom assertthat has_extension
#' @importFrom assertthat is.string
#' @importFrom assertthat has_attr
#' @importFrom assertthat "on_failure<-"
#' @importFrom assertthat not_empty
#' @importFrom assertthat noNA
NULL

## NCBI BLAST defline database tags
.tags  <- c("lcl","gb","emb","pir","sp","ref","gnl", "gi","dbj","prf","pdb","pat","bbs")

is.empty <- function(x) {
  is.null(x) || length(x) == 0L || (length(x) == 1L && !nzchar(x))
}
on_failure(is.empty) <- function(call, env) {
  paste0(deparse(call$x), " is not empty.")
}

"%||%" <- function(a, b) {
  if (is.empty(a)) force(b) else a
}

"%|na|%" <- function(a, b) {
  if (is.null(a) || all(is.na(a))) force(b) else a
}

## Vectorized default operators
"%|%" <- function(a, b) ifelse(nzchar(a), a, b)

"%|NA|%" <- function(a, b) ifelse(is.na(a), b, a)

"%ni%" <- Negate(`%in%`)

compact <- function(x) {
  x[!vapply(x, is.null, FALSE, USE.NAMES=FALSE)]
}

compactChar <- function(x) {
  x[vapply(x, nzchar, FALSE, USE.NAMES=FALSE)]
}

Call <- function(fn, ...) {
  fn <- match.fun(fn)
  fn(...)
}

Partial <- function(fn, ..., .env = parent.frame()) {
  assert_that(is.function(fn))
  fcall <- substitute(fn(...))
  if (!is.primitive(fn))
    fcall <- match.call(fn, fcall)  
  fcall[[length(fcall) + 1]] <- quote(...)
  args <- list("..." = quote(expr = ))
  eval(call("function", as.pairlist(args), fcall), .env)
}

Compose <- function (...) {
  fns <- lapply(compact(list(...)), match.fun)
  len <- length(fns)
  function (...) {
    res <- Call(fns[[len]], ...)
    for (fn in rev(fns[-len]))
      res <- fn(res)
    res
  }
}

merge_list <- function(x, y) {
  if (length(x) == 0) return(y)
  if (length(y) == 0) return(x) 
  i <- is.na(match(names(y), names(x)))
  if (any(i)) {
    x[names(y)[which(i)]] = y[which(i)]
  }
  x
}

trim <- function(x, trim = '\\s+') {
  assert_that(is.vector(x))
  gsub(paste0("^", trim, "|", trim, "$"), '', x)
}

dup <- function (x, n) {
  assert_that(is.string(x))
  if (any(n < 0)) n[n < 0] <- 0
  vapply(.mapply(rep.int, list(rep.int(x, length(n)), n), NULL), paste0, collapse="", "")
}

blanks <- Partial(dup, x = " ")

pad <- function (x, n = 10, where = 'left', pad = ' ') {
  assert_that(length(n) == 1, length(pad) == 1)
  x <- as.character(x)
  where <- match.arg(where, c("left", "right", "both"))
  needed <- pmax(0, n - nchar(x))
  left <- switch(where, left = needed, right = 0, both = floor(needed/2))
  right <- switch(where, left = 0, right = needed, both = ceiling(needed/2))
  lengths <- unique(c(left, right))
  padding <- dup(pad, lengths)
  paste0(padding[match(left, lengths)], x, padding[match(right, lengths)])
}

ellipsize <- function(obj, width = getOption("width"), ellipsis = " ...") {
  str <- encodeString(obj)
  ifelse(nchar(str) > width - 1,
         paste0(substring(str, 1, width - nchar(ellipsis)), ellipsis),
         str)
}

strsplitN <- function (x, split, n, from = "start", collapse = split, ...) {
  assert_that(is.vector(x))
  from <- match.arg(from, c("start", "end"))
  xs <- strsplit(x, split, ...)
  end <- vapply(xs, length, integer(1))
  
  if (from == "end") {
    end <- end + 1L
    n <- lapply(end, `-`, n)
    n <- .mapply(`[<-`, list(x=n, i=lapply(n, `<`, 0), value=0L), NULL)
  } else {
    n <- lapply(rep(0, length(xs)), `+`, n)
    n <- .mapply(`[<-`, list(x=n, i=Map(`>`, n, end), value=end), NULL)
  }  
  n <- lapply(n, Compose("sort", "unique"))
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse), list(x = xs, n = n), NULL))
}

#' Format paragraphs
#' 
#' Similar to \code{\link{strwrap}} but returns a single string with
#' linefeeds inserted
#' 
#' @param s a character vector or a list of character vectors
#' @param width a positive integer giving the column for inserting
#' linefeeds
#' @param indent an integer giving the indentation of the first line of
#' the paragraph; negative values of \code{indent} are allowed and reduce
#' the width for the first line by that value.
#' @param offset a non-negative integer giving the indentation of all
#' but the first line
#' @param split regular expression used for splitting. Defaults to
#' a whitespace character.
#' @param FORCE Words are force split if the available width is too small.
#' @param FULL_FORCE Lines are split exactly at the specified width
#' irrespective of whether there is whitespace or not.
#' 
#' @return a character vector
#' @keywords internal
linebreak <- function(s, width = getOption("width") - 2,
                      indent = 0, offset = 0, split = " ",
                      FORCE = FALSE, FULL_FORCE = FALSE)
{
  if (!is.character(s)) 
    s <- as.character(s)
  
  if (length(s) == 0L)
    return("")
  
  .first_iteration <- TRUE
  # set indent string to "" if a negative value is given
  # this lets us shrink the available width for the first line by that value
  indent_string <- blanks(indent)
  offset_string <- paste0("\n", blanks(offset))
  
  ans <- Map(function(s, width, offset, indent,
                      indent_string, split, FORCE,
                      FULL_FORCE)
  {
    # remove leading and trailing blanks
    # convert newlines, tabs, spaces to " "
    # find first position where 'split' applies
    if (!FULL_FORCE) {
      s <- gsub("\\s+", " ", trim(s), perl=TRUE)
    }
    fws <- regexpr(split, s, perl=TRUE)
    
    if (.first_iteration)
      string_width <- indent + nchar(s)
    else
      string_width <- offset + nchar(s)
    
    if (string_width > width)
    {
      # if not everything fits on one line
      .first_iteration <- FALSE
      if(FULL_FORCE ||
           ((fws == -1 || fws >= (width - string_width)) && FORCE))
      {
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
      else if ((fws == -1 || fws >= (width - offset + indent)) && !FORCE)
      {
        # if no whitespace or first word too long and NO force break
        # stop right here
        stop("Can't break in the middle of a word. Use the force!")
      }
      else
      {
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
    } else {
      # if everything fits on one line go with the string + indent
      paste0(indent_string, s)
    }
  }, s, width, offset, abs(indent), indent_string, split,
             FORCE, FULL_FORCE, USE.NAMES=FALSE)
  unlist(ans)
}

xvalue <- function(doc, path, as = 'character', default = NA_character_, fun = NULL, ...) {
  fun <- Compose(fun, xmlValue)
  res <- unlist(xpathApply(doc, path, fun, ...)) %||% default
  set_mode(res, as)
}

xsize <- function(doc, path, ...) {
  length(xpathApply(doc, path, ...))
}

set_mode <- function(x, as) {
  AS <- match.fun(paste0('as.', as))
  if (!is.null(x)) AS(x) else x
}

#' @keywords internal
#' @export
subl <- function(x, ...) {
  assert_that(has_command('subl'))
  if (tryCatch(is.readable(x), assertError = function (e) FALSE)) {
    SysCall('subl', stdin=x, redirection=FALSE, ...)
  }
  else {
    tmp <- tempfile()
    write(x, file=tmp)
    SysCall('subl', stdin=tmp, redirection=FALSE, ...)
  }
}

are_null <- function (x) {
  vapply(x, is.null, FALSE, USE.NAMES=FALSE)
}

are_true <- function (x) {
  vapply(x, isTRUE, FALSE, USE.NAMES=FALSE)
}

are_false <- function (x) {
  vapply(x, function (x) identical(x, FALSE), FALSE, USE.NAMES=FALSE)
}

#' Test if an external executable is available
#' 
#' Uses \code{\link{Sys.which}} internally, so it should work
#' on Windows and Unix.alikes.
#' 
#' @param cmd The exececutable to test for.
#' @param msg Additional message if the test fails.
#' @keywords internal
has_command <- function (cmd, msg = "") {
  assert_that(is.string(cmd))
  unname(Sys.which(cmd) != "")
}
on_failure(has_command) <- function(call, env) {
  paste0("Dependency ", sQuote(eval(call$cmd, env)), " is not installed\n",
         eval(call$msg, env))
}

#' Wrapper for system commands
#' 
#' @param exec The system command to be invoked.
#' @param ... Arguments passed on to the \code{system} command as name-value or 
#' name=\code{TRUE} pairs.
#' @param args Named list of arguments passed on to the \code{system} command.
#' Is merged with \code{...}.
#' @param stdin Input.
#' @param stdout Output.
#' @param redirection Redirection.
#' @param style One of \sQuote{unix} or \sQuote{gnu}.
#' @param sep Seperator of option and option argument.
#' @param show_cmd Have a look what the final command looks like.
#' @param intern Passed on to \code{\link{system}}'s \code{intern} argument.
#' @param input Passed on to \code{\link{system}}'s \code{input} argument.
#' @keywords internal
SysCall <- function (exec, ..., args = list(), stdin = NULL, stdout = NULL,
                     redirection = TRUE, style = c("unix", "gnu"), sep = " ",
                     show_cmd = FALSE, intern = FALSE, input = NULL) {  
  assert_that(has_command(exec))
  args <- merge_list(list(...), args)
  style <- match.arg(style)
  if (is.null(stdin)) {
    stdin <- ""
  } else if (!is.null(stdin) && redirection) {
    stdin <- paste("<", stdin)
  }
  if (is.null(stdout)) {
    stdout <- ""
  } else {
    stdout <- paste(">", stdout)
  }
  args[are_true(args)] <- ""
  args[are_false(args) | are_null(args)] <- NULL
  args <- switch(style,
                 unix=paste0(trim(sprintf("-%s%s%s", names(args), sep, args)), collapse=" "),
                 gnu=paste0(trim(sprintf("--%s%s%s", names(args), sep, args)), collapse=" "))
  
  if (show_cmd)
    print(trim(paste(exec, args, stdin, stdout)))
  else
    system(trim(paste(exec, args, stdin, stdout)), intern = intern, input = input)
}

