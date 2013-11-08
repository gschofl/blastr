#' @include sqlite-utils.r
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
  function(x, id) {
    args <- list(...)
    stmts <- trim(paste("SELECT", SELECT, "FROM", FROM,
                        if (!is.null(args$WHERE)) {
                          paste("WHERE", args$WHERE, "=", id)
                        },
                        if (!is.null(args$VAL) && !is.null(args$FUN)) {
                          paste("AND", args$VAL, "= (SELECT", args$FUN,
                                "(", args$VAL, ") FROM", FROM, "WHERE", args$WHERE, "=", id, ")")
                        }))
    AS <- match.fun(paste0('as.', as))
    lapply(stmts, function(stmt) {
      AS(db_query(x, stmt, 1L) %||% NA_character_)
    })
  }
}

getterFromToRange <- function(x, id, type='query', max=FALSE) {
  if (max) {
    if (type=='query') {
      pos <- db_query(x,paste('SELECT hit_id, query_frame, query_from, query_to from hsp 
                              WHERE query_id=',id, 'AND bit_score = (SELECT
                              MAX(bit_score) FROM hsp WHERE query_id =', id, ')'))
    } else {
      pos <- db_query(x,paste('SELECT hit_id, hit_frame, hit_from, hit_to from hsp 
                              WHERE query_id=',id, 'AND bit_score = (SELECT
                              MAX(bit_score) FROM hsp WHERE query_id =', id, ')'))
    }
  } else {
    if (type=='query') {
      pos <- db_query(x,paste('SELECT hit_id, query_frame, query_from, query_to from hsp 
                                WHERE query_id =',id))
    } else {
      pos <- db_query(x,paste('SELECT hit_id, hit_frame, hit_from, hit_to from hsp 
                                WHERE query_id =',id))
    }
  }
  pos
}
