library(rmisc)
library(RSQLite)

con <- blastReportDB("/home/psehnert/daten/SPICEIII/miseq/sample64/blast.test.db")
con

#query Table
getQueryDef(con,9)
getQueryLen(con,9)

# hit table
getHitID(con,9)
getHitNum(con,9)
getHitLen(con,9)
getAccession(con,9)
getGeneID(con,9)
getHitDef(con,9)
getHit(con,9)

getHspHitID(con,20,max=T)
getHspID(con,20)
getHspNum(con,20)
getBitscore(con,9)
getScore(con,9)
getEvalue(con,9)
getQueryFrom(con,9)
getQueryTo(con,9)
getQueryFrame(con,9)
getHitFrom(con,9,max=T)
getHitTo(con,9)
getHitFrame(con,9)
getIdentity(con,9)
getPositive(con,9)
getGaps(con,9)
getAlignLen(con,9)
getQuerySeq(con,9)
getHitSeq(con,9)
getMatch(con,9)
getPercIdentity(con,9)
getHitCoverage(con,20)

as.list( unlist(getIdentity(x, id))/unlist(getAlignLen(x, id))
         
db_query(con,"SELECT identity,align_len from hsp where query_id=9")

db_query(con,"SELECT CAST(identity AS FLOAT)/CAST(align_len AS FLOAT) as PERC
              FROM (SELECT identity,align_len 
                    FROM hsp 
                    WHERE query_id=9)", 1)
        
         
         
getIdentity(con,9)
getAlignLen(con,9)
