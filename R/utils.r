#' @importFrom assertthat is.string
#' @importFrom ape read.dna
#' @importFrom ape write.dna
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readDNAMultipleAlignment
#' @import methods
NULL

setOldClass("DNAbin")

## NCBI BLAST defline database tags
.tags  <- c("lcl","gb","emb","pir","sp","ref","gnl",
            "gi","dbj","prf","pdb","pat","bbs")

## vectorised %||%
"%|%" <- function (a, b) {
  ifelse(is_empty(a), b, a)
}

## vectorised %|na|%
"%|NA|%" <- function (a, b) {
  ifelse(is.na(a), b, a)
}

listclassConstructor <- function (listClass, elemClass) {
  assert_that(is.string(listClass), is.string(elemClass))
  function (..., query_env) {
    listData <- list(...)
    if (length(listData) == 0L) {
      new(listClass, list(new(elemClass)))
    } else {
      if (length(listData) == 1L && is.list(listData[[1L]])) 
        listData <- listData[[1L]]
      if (!all(vapply(listData, is, elemClass, FUN.VALUE=logical(1L)))) 
        stop("All elements in '...' must be '", elemClass,"' objects")
      if (!missing(query_env)) {
        new(listClass, .Data = listData, query_env = query_env)
      } else {
        new(listClass, .Data = listData)
      }
      
    }
  }
}

listclassValidator <- function (listClass, elemClass) {
  assert_that(is.string(listClass), is.string(elemClass))
  function (object) {
    errors <- character()
    elem_of_class <- vapply(S3Part(object, strictS3=TRUE), is, elemClass, FUN.VALUE=logical(1L))
    if (!all(elem_of_class)) {
      msg <- paste0("All elements in a '", listClass ,"' must be of class '",
                    elemClass, "'.")
      errors <- c(errors, msg)
    }
    
    if (length(errors) == 0L) TRUE else errors
  }
}

getterConstructor <- function(SELECT, FROM, ..., as = 'character') {
  function (x, id, verbose = FALSE) {
    args <- list(...)
    stmts <- trim(paste("SELECT", SELECT, "FROM", FROM,
                        if (!is.null(args$WHERE)) {
                          paste("WHERE", args$WHERE, "=", id)
                        },
                        if (!is.null(args$VAL) && !is.null(args$FUN)) {
                          paste("AND", args$VAL, "= (SELECT", args$FUN,
                                "(", args$VAL, ") FROM", FROM, "WHERE", args$WHERE, "=", id, ")")
                        }))

    if (verbose) cat(stmts)
    AS <- match.fun(paste0('as.', as))
    lapply(stmts, function(stmt) {
      AS( db_query(x, stmt, 1L) %||% NA_character_ )
    })
  }
}

# advancedGetterConstructor <- function(SELECT,FROM,WHERE,SQL_FUNCTION,VALUE) {
#   function(x,id) {
#     lapply(id,function(id) {
#       SQL <- paste("SELECT", SELECT, "FROM", FROM, "WHERE", WHERE, "=",id,
#                    "AND", VALUE, "= (SELECT",SQL_FUNCTION,"(",VALUE,") FROM", 
#                     FROM, "WHERE",WHERE, "=",id,")")
#       db_query(x,SQL,1L) %||% NA_character_
#     })
#   }
# }

ellipsize <- function(obj, width = getOption("width"), ellipsis = " ...") {
  str <- encodeString(obj)
  ifelse(nchar(str) > width - 1,
         paste0(substring(str, 1, width - nchar(ellipsis)), ellipsis),
         str)
}

setAs("XStringSet", "DNAbin", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeXStringSet(from, filepath=tmp, format="fasta")
  to <- read.dna(tmp, format="fasta")
  to
})


setAs("DNAbin", "DNAStringSet", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  write.dna(from, file=tmp, format="fasta")
  to <- readDNAStringSet(tmp, "fasta")
  to
})


setAs("DNAbin", "DNAMultipleAlignment", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  write.dna(from, file=tmp, format="fasta")
  to <- readDNAMultipleAlignment(tmp, "fasta")
  to
})


setAs("XStringSet", "DNAMultipleAlignment", function (from) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeXStringSet(from, filepath=tmp, format="fasta")
  to <- readDNAMultipleAlignment(tmp, "fasta")
  to
})

