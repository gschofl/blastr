
# Blast-Getters -------------------------------------------------------------

#' @importFrom IRanges IRanges
#' @importFrom IRanges RangedData
NULL

# query ============================
#
#' Getter methods for blast records
#' 
#' @param x A \code{\linkS4class{blastRecord}} object.
#' @param ... Additional arguments
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("query", function (x, ...) standardGeneric("query"))


#' @export
setMethod("query", "blastReport", function (x) x@query)


# index ============================
#
#' Getter methods for blast records
#' 
#' @param x A \code{\linkS4class{blastRecord}} object.
#' @param attributes Include definition as a name attribute.
#' 
#' @rdname blastReport-method
#' @export
setMethod("index", "blastReport",
          function (x, attributes = FALSE) {
            if (is.na(suppressWarnings(as.numeric(x@query$identifier)))) {
              message("No index")
              return(invisible(FALSE))
            }
            if (attributes) {
              setNames(as.numeric(x@query$identifier), nm=x@query$def)
            } else {
              as.numeric(x@query$identifier)
            }
          })


# hits ============================
#
#' Getter methods for blast records
#' 
#' @param x A \code{\linkS4class{blastRecord}}  object.
#' @param ... Additional arguments
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("hits", function (x, ...) standardGeneric("hits"))

#' @export
setMethod("hits", "blastReport", function (x) x@hits)


# getId ===========================
#
#' Getter methods for blast records
#' 
#' @usage getId(x, db="gi")
#' 
#' @param x A \code{\linkS4class{blastRecord}} or
#' \code{\linkS4class{hit}} object.
#' @param db Database tag (e.g.: 'gi', 'gb', 'emb', 'ref', 'gnl', ...)
#'
#' @return Accession numbers as specified by \code{db} for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getId", function (x, db="gi", ...) standardGeneric("getId"))


#' @export
setMethod("getId", "blastReport", function (x, db, ...) {
  lapply(x@hits, function (x) getId(x, db)) 
})


#' @export
setMethod("getId", "hit", function (x, db, ...) {
  id <- lapply(x@id, "[[", db)
  id[vapply(id, is.null, logical(1))] <- NA_character_
  unlist(id)
})


## getDesc =========================
##
#' Getter methods for blast records
#' 
#' @usage getDesc(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Description for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getDesc", function (x) standardGeneric("getDesc"))


#' @export
setMethod("getDesc", "blastReport", function (x) {
  lapply(x@hits, function (x) getDesc(x)) 
})


#' @export
setMethod("getDesc", "hit", function (x) {
  unlist(x@desc)
})


## getAccn =========================
##
#' Getter methods for blast records
#' 
#' @usage getAccn(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Accession number for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getAccn", function (x) standardGeneric("getAccn"))


#' @export
setMethod("getAccn", "blastReport", function (x) {
  lapply(x@hits, function (x) getAccn(x)) 
})


#' @export
setMethod("getAccn", "hit", function (x) unlist(x@accn) )


## getBitscore =====================
##
#' Getter methods for blast records
#' 
#' @usage getBitscore(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Bitscores for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getBitscore", function (x) standardGeneric("getBitscore"))


#' @export
setMethod("getBitscore", "blastReport", function (x) {
  lapply(x@hits, function (x) getBitscore(x)) 
})


#' @export
setMethod("getBitscore", "hit", function (x) unlist(x@hsp@bit_score) )


## getEvalue =======================
##
#' Getter methods for blast records
#' 
#' @usage getEvalue(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Evalue for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getEvalue", function (x) standardGeneric("getEvalue"))


#' @export
setMethod("getEvalue", "blastReport", function (x) {
  lapply(x@hits, function (x) getEvalue(x)) 
})


#' @export
setMethod("getEvalue", "hit", function (x) unlist(x@hsp@evalue) )


## getQueryRange ===================
##
#' Getter methods for blast records
#' 
#' @usage getQueryRange(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Start and end position for each query.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getQueryRange", function (x) standardGeneric("getQueryRange"))


#' @export
setMethod("getQueryRange", "hit", function (x) {
  frame <- x@hsp@query_frame
  start <- as.integer(ifelse(frame >= 0, x@hsp@query_from, x@hsp@query_to))
  end <- as.integer(ifelse(frame >= 0, x@hsp@query_to, x@hsp@query_from))
  RangedData(IRanges(start, end), frame = frame)
}) 


#' @export
setMethod("getQueryRange", "blastReport", function (x) {
  lapply(x@hits, getQueryRange) 
})


## getHitRange =====================
##
#' Getter methods for blast records
#' 
#' @usage getHitRange(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Start and end position for each hit.
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getHitRange", function (x) standardGeneric("getHitRange"))


#' @export
setMethod("getHitRange", "hit", function (x) {
  frame <- x@hsp@hit_frame
  start <- as.integer(ifelse(frame >= 0, x@hsp@hit_from, x@hsp@hit_to))
  end <- as.integer(ifelse(frame >= 0, x@hsp@hit_to, x@hsp@hit_from))
  RangedData(IRanges(start, end), frame = frame)
})



#' @export
setMethod("getHitRange", "blastReport", function (x) {
  lapply(x@hits, getHitRange) 
})


## getHitFrame =====================
##
#' Getter methods for blast records
#' 
#' @usage getHitFrame(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Frame for each hit (1 or -1)
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getHitFrame", function (x) standardGeneric("getHitFrame"))


#' @export
setMethod("getHitFrame", "blastReport", function (x) {
  lapply(x@hits, function (x) getHitFrame(x)) 
})


#' @export
setMethod("getHitFrame", "hit", function (x) x@hsp@hit_frame)


## getQuerySeq =====================
##
#' Getter methods for blast records
#' 
#' @usage getQuerySeq(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Query sequence(s) as a \code{\link[Biostrings]{DNAStringSet}}
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getQuerySeq", function (x) standardGeneric("getQuerySeq"))


#' @export
setMethod("getQuerySeq", "blastReport", function (x) {
  lapply(x@hits, function (x) getQuerySeq(x)) 
})


#' @export
setMethod("getQuerySeq", "hit", function (x) x@hsp@qseq)


# getHitSeq =======================
#
#' Getter methods for blast records
#' 
#' @usage getHitSeq(x)
#' 
#' @param x A \code{\link{blastReport-class}} or \code{\link{hit-class}} object.
#' 
#' @return Hit sequence(s) as a \code{\link[Biostrings]{DNAStringSet}}
#' 
#' @rdname blastReport-method
#' @docType methods
#' @export
setGeneric("getHitSeq", function (x) standardGeneric("getHitSeq"))


#' @export
setMethod("getHitSeq", "blastReport", function (x) {
  lapply(x@hits, function (x) getHitSeq(x)) 
})


#' @export
setMethod("getHitSeq", "hit", function (x) x@hsp@hseq)

# BlastDB-Getters -------------------------------------------------------------

## Query Table

setGeneric("getQueryDef", function (con, query_id) standardGeneric("getQueryDef"))


#' @export
setMethod("getQueryDef", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT query_def 
                  FROM query 
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getQueryLen", function (con, query_id) standardGeneric("getQueryLen"))

#' @export
setMethod("getQueryLen", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT query_len 
                  FROM query 
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

## Hit Table

setGeneric("getHitID", function (con, query_id) standardGeneric("getHitID"))

#' @export
setMethod("getHitID", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_id
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitNum", function (con, query_id) standardGeneric("getHitNum"))

#' @export
setMethod("getHitNum", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_num
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getGeneID", function (con, query_id) standardGeneric("getGeneID"))

#' @export
setMethod("getGeneID", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT gene_id
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getAccn", function (con, query_id) standardGeneric("getAccn"))

#' @export
setMethod("getAccn", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT accession
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitDef", function (con, query_id) standardGeneric("getHitDef"))

#' @export
setMethod("getHitDef", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT definition
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitLen", function (con, query_id) standardGeneric("getHitLen"))

#' @export
setMethod("getHitLen", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT length
                  FROM hit
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

## Hsp Table

setGeneric("getHspHitID", function (con, query_id) standardGeneric("getHspHitID"))

#' @export
setMethod("getHspHitID", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_id
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHspID", function (con, query_id) standardGeneric("getHspID"))

#' @export
setMethod("getHspID", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hsp_id
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHspNum", function (con, query_id) standardGeneric("getHspNum"))

#' @export
setMethod("getHspNum", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_num
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getBitscore", function (con, query_id) standardGeneric("getBitscore"))

#' @export
setMethod("getBitscore", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT bit_score
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getScore", function (con, query_id) standardGeneric("getScore"))

#' @export
setMethod("getScore", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT score
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getEvalue", function (con, query_id) standardGeneric("getEvalue"))

#' @export
setMethod("getEvalue", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT evalue
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getQueryStart", function (con, query_id) standardGeneric("getQueryStart"))

#' @export
setMethod("getQueryStart", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT query_from
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getQueryEnd", function (con, query_id) standardGeneric("getQueryEnd"))

#' @export
setMethod("getQueryEnd", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT query_to
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitStart", function (con, query_id) standardGeneric("getHitStart"))

#' @export
setMethod("getHitStart", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_from
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitEnd", function (con, query_id) standardGeneric("getHitEnd"))

#' @export
setMethod("getHitEnd", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_to
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})



setGeneric("getQueryFrame", function (con, query_id) standardGeneric("getQueryFrame"))

#' @export
setMethod("getQueryFrame", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT query_frame
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getHitFrame", function (con, query_id) standardGeneric("getHitFrame"))

#' @export
setMethod("getHitFrame", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hit_frame
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getIdentity", function (con, query_id) standardGeneric("getIdentity"))

#' @export
setMethod("getIdentity", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT identity
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

setGeneric("getPositive", function (con, query_id) standardGeneric("getPositive"))

#' @export
setMethod("getPositive", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT positive
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})


setGeneric("getGaps", function (con, query_id) standardGeneric("getGaps"))

#' @export
setMethod("getGaps", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT gaps
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})


setGeneric("getAlignLen", function (con, query_id) standardGeneric("getAlignLen"))

#' @export
setMethod("getAlignLen", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT align_len
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})


setGeneric("getQSeq", function (con, query_id) standardGeneric("getQSeq"))

#' @export
setMethod("getQSeq", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT qseq
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})


setGeneric("getHSeq", function (con, query_id) standardGeneric("getHSeq"))

#' @export
setMethod("getHSeq", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT hseq
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})


setGeneric("getMidline", function (con, query_id) standardGeneric("getMidline"))

#' @export
setMethod("getMidline", "blastReportDB", function (con,query_id) {
  lapply(query_id, function(x) {
    sql <- paste0("SELECT midline
                  FROM hsp
                  WHERE query_id = ", 
                  x)
    db_query(con, sql, 1L)
  })
})

