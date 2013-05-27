require(XML)
require(blastr)
require(RSQLite)


###Aufgaben

# klasse erstellen blastdb mit con
setClass('blastReportDB', contains='SQLiteConnection')

setValidity('blastReportDB', function (object) {
  if (!all(c("hit","hsp","query") %in% dbListTables(object)))
    return("Table missing from 'blastReportDB'")
  if (!all(c("query_id","query_def","query_len") %in% dbListFields(object, "query")))
    return("Field missing from table 'query'")
  if (!all(c("query_id","hit_id","hit_num","gene_id","accession",
             "definition","length") %in% dbListFields(object, "hit")))
    return("Field missing from table 'query'")
  if (!all(c("query_id","hit_id","hsp_id","hsp_num","bit_score",
             "score","evalue","query_from","query_to","hit_from",
             "hit_to","query_frame","hit_frame","identity","positive",
             "gaps","align_len","qseq","hseq","midline") 
           %in% dbListFields(object, "hsp")))   
    return("Field missing from table 'query'")
  TRUE
})


# getter für jeden slot der table

# methoden für coverage query und hit
# metode für perc_ident

# methode für connect (blastReportDB)
#  - connection erstellen
# methode für disconect



# connect do db
require(rmisc)




con <- blastReportDB("/home/psehnert/daten/SPICEIII/miseq/sample64/blast.test.db")

con
blastReportDB <- function( dbname ) {
  con <- db_connect(dbName=dbname, message="")
  new('blastReportDB', con)
}

# show methode ausbauen für object beschreibung
setMethod("show",
          signature(object = "blastReportDB"),
          function (object) 
          {
            n <- db_count(object, "query")
            showme <- sprintf('%s object with %s query rows',
                              sQuote(class(object)), n)
            cat(showme, sep="\n")
          }
)

setMethod("show", "blastReport",
         function (object) {
           def_pos <- which(names(object@query) ==  "def")
           def_line <- deparseDeflines(list(object@query[1:def_pos - 1]),
                                       object@query[3])
           query <- linebreak(sprintf("%s (%s letters)",
                                      def_line, object@query[["len"]]),
                              offset = 10)     
           cat(sprintf("Query:    %s\nProgram:  %s\nDatabase: %s\n\n",
                       query, sQuote(object@version), sQuote(object@db)))
           cat(paste0("Accession      ",
                      format("Description", width=ceiling(getOption("width")*0.5)),
                      " bit score", "  evalue\n"))
           invisible(
             lapply(object@hits, function (x) {
               cat(paste(format(x@accn, width=14),
                         format(
                           strtrim(deparseDeflines(x@id, x@desc)[[1]],
                                   width=floor(getOption("width")*0.5)),
                           width=ceiling(getOption("width")*0.5)),
                         format(x@hsp@bit_score, digits=4, width=9),
                         format(x@hsp@evalue, scientific=TRUE, digits=2, width=8),
                         "\n"))
             }))
           
           return(invisible(NULL))
         })






con <- blastReportDB("/home/psehnert/daten/SPICEIII/miseq/sample64/blast.test.db")
con




