require(XML)
require(blastr)
require(RSQLite)

# connect do db

con <- dbConnect(SQLite(),dbname="/home/psehnert/daten/SPICEIII/miseq/sample64/blast.test.db")
res <- dbGetQuery(con, "SELECT * FROM query")
res <- fetch(dbSendQuery(con,"SELECT * FROM query"),n=-1)

