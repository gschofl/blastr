require(XML)
require(blastr)
require(RSQLite)


###Aufgaben

# klasse erstellen blastdb mit con

# getter für jeden slot der table

# methoden für coverage query und hit
# metode für perc_ident

# methode für connect (blastReportDB)
#  - connection erstellen
#  - validität (alle felder alle tables) checken
# methode für disconect



# connect do db
require(rmisc)

setClass('blastReportDB', contains='SQLiteConnection')

setValidity('blastReportDB', function (object) {
  if (!all(c("hit","hsp","query") %in% dbListTables(object)))
    return("Table missing from 'blastReportDB'")
  if (!all(c("query_id","query_def","query_len","bla") %in% dbListFields(object, "query")))
    return("Field missing from table 'query'")
  
  
  TRUE
})

blastReportDB <- function( dbname ) {
  con <- db_connect(dbName=dbname, message="")
  new('blastReportDB', con)
}

# einfache art und weise sich schon skellte für getter zu erstellen
method.skeleton("show", signature="blastReportDB")
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







con <- blastReportDB("/home/psehnert/daten/SPICEIII/miseq/sample64/blast.test.db")
con


# query table
query <- fetch(dbSendQuery(con,"SELECT * FROM query"),n=-1)
query_id <- fetch(dbSendQuery(con,"SELECT query_id FROM query WHERE query_id = 1"),n=-1)
query_def <- fetch(dbSendQuery(con,"SELECT query_def FROM query WHERE query_id = 1"),n=-1)
query_len <- fetch(dbSendQuery(con,"SELECT query_len FROM query WHERE query_id = 1"),n=-1)

# hit table
hit <- fetch(dbSendQuery(con,"SELECT *  FROM hit"),n=-1)
query_id <-  fetch(dbSendQuery(con,"SELECT query_id FROM hit WHERE query_id = 9"),n=-1)
hit_id	<- fetch(dbSendQuery(con,"SELECT hit_id FROM hit WHERE query_id = 9"),n=-1)
hit_num	<- fetch(dbSendQuery(con,"SELECT hit_num FROM hit WHERE query_id = 9"),n=-1)
gene_id	<- fetch(dbSendQuery(con,"SELECT gene_id  FROM hit WHERE query_id = 9"),n=-1)
accession	<- fetch(dbSendQuery(con,"SELECT accession  FROM hit WHERE query_id = 9"),n=-1)
definition <- fetch(dbSendQuery(con,"SELECT definition FROM hit WHERE query_id = 9"),n=-1)	
length <- fetch(dbSendQuery(con,"SELECT length  FROM hit WHERE query_id = 9"),n=-1)

# hsp table
hsp <- fetch(dbSendQuery(con,"SELECT *  FROM hsp"),n=-1)