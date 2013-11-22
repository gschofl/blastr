#' @include sqlite-db.r
NULL


listclassConstructor <- function(listClass, elemClass) {
  assert_that(is.string(listClass), is.string(elemClass))
  function(..., query_env) {
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


listclassValidator <- function(listClass, elemClass) {
  assert_that(is.string(listClass), is.string(elemClass))
  function(object) {
    errors <- character()
    elem_of_class <- vapply(S3Part(object, strictS3=TRUE), is, elemClass, FUN.VALUE=logical(1L))
    if (!all(elem_of_class)) {
      msg <- paste0("All elements in a '", listClass ,"' must be of class '", elemClass, "'.")
      errors <- c(errors, msg)
    }
    if (length(errors) == 0L) TRUE else errors
  }
}

getterConstructor <- function(SELECT, FROM, ..., as = 'character') {
  ARGS <- list(...)
  function(x, id, ...) {
    AS <- match.fun(paste0('as.', as))
    con <- conn(x)
    log <- list(...)$log
    WHERE <- AND <- GROUPBY <- ''
    if (!is.null(ARGS$WHERE) && !missing(id)) {
      WHERE <- paste0(" where ", ARGS$WHERE, " = ", id) 
    }
    if (!is.null(ARGS$GROUP_BY)) {
      GROUPBY <- paste0(" group by ", ARGS$GROUP_BY)
    }
    if (!is.null(ARGS$VAL) && !is.null(ARGS$FUN)) {
      AND <-  paste0(if (all(nzchar(WHERE))) " and " else " where ",
                         ARGS$VAL, " in (select ", ARGS$FUN,
                         "(", ARGS$VAL, ") from ", FROM, WHERE, GROUPBY, ")")
    }
    stmts <- trim(paste0("select ", SELECT, " from ", FROM, WHERE, AND, GROUPBY))
    lapply(stmts, function(stmt) {
      AS(.db_query(con, stmt, 1L, log=log) %||% NA)
    })
  }
}


simpleGetter <- function(SELECT, FROM, ..., as = 'character') {
  ARGS <- list(...)
  function(x, id, ...) {
    AS <- match.fun(paste0('as.', as))
    log <- list(...)$log
    WHERE <- ''
    if (!is.null(ARGS$WHERE) && !missing(id)) {
      WHERE <- paste0(" where ", ARGS$WHERE, " in (", paste0(id, collapse=","), ")") 
    }
    stmt <- trim(paste0("select ", SELECT, " from ", FROM, WHERE))
    AS(.db_query(conn(x), stmt, 1L, log=log) %||% NA)
  }
}


.rangeDB <- function(con, id, type, width=FALSE, max=FALSE, ...) {
  log <- list(...)$log
  stmts <- paste0("select hit_id as id, min(", type, "_from, ", type, "_to) as start,",
                  " abs(", type, "_to - ", type, "_from + 1) as width",
                  " from hsp where query_id = ", id,
                  if (max) 
                    paste0(' and bit_score = (select max(bit_score)',
                           ' from hsp where query_id = ', id, ')')
                  else '')
  pos <- lapply(stmts, .db_query, con=con, j=NA, log=log)
  lapply(pos, function(pos) {
    lapply(unname(split.data.frame(pos, as.factor(pos$id))), function(p) {
      r <- IRanges(start=p$start, width=p$width)
      if (width) width(reduce(r)) else r
    })
  })
}


