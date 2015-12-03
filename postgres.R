# install.packages("devtools")
# devtools::install_github("RcppCore/Rcpp")
# devtools::install_github("rstats-db/DBI")
# devtools::install_github("rstats-db/RPostgres")

library(DBI)
# Connect to our postgres database with a readonly user
con <- dbConnect(RPostgres::Postgres(),dbname = 'toshi', 
                 host = 'toshi.cn6zzwcfsto5.us-east-1.rds.amazonaws.com',
                 port = 5432,
                 user = 'readonly',
                 password = 'password')
dbListTables(con)

#construct a query
res <- dbSendQuery(con, "SELECT * FROM blocks WHERE height = 0")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}

#fetch the results of the query
dbFetch(res)
dbClearResult(res)


#let's count how many transactions are in our db
res <- dbSendQuery(con, "select count(*) from transactions")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
dbFetch(res)
dbClearResult(res)


#let's count how many blocks are in our db
res <- dbSendQuery(con, "select count(*) from blocks")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
dbFetch(res)
dbClearResult(res)


#let's see the timestamp of the most recent block in our database
res <- dbSendQuery(con, "select time, height from blocks order by height desc limit 1")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
time <- dbFetch(res)
print(time$height[1])
print(as.POSIXct(time$time[1], origin="1970-01-01"))
dbClearResult(res)


#let's calculate some block stats
res <- dbSendQuery(con, "select height, time, size, transactions_count from blocks")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
blockstats <- dbFetch(res)
dbClearResult(res)

# I'm sure this is horrible code...  
blockstats$elapsed_time <-  blockstats$time - c(NA,blockstats$time[-length(blockstats$time)])
## cleaner but slower: c(NA, diff(blockstats$time))
# for some reason this results in some negative values for elapsed time. not sure why that might be
# maybe some of the computers that discovered blocks had clocks that were off? by up to 2 hours?
# might be a problem for only early in the chain, though

blockstats$txpermin <- blockstats$transactions_count*60/blockstats$elapsed_time
blockstats$bytespermin <- blockstats$size*60/blockstats$elapsed_time

library(Hmisc)
describe(blockstats)


#let's calculate some block stats for the last 1000 blocks in our db
res <- dbSendQuery(con, "select height, time, size, transactions_count from blocks order by height desc limit 1000")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
blockstats1k <- dbFetch(res)
dbClearResult(res)
blockstats1k <-blockstats1k[order(blockstats1k$height),] # re-sorting because we selected blocks in reverse order
blockstats1k$elapsed_time <- blockstats1k$time - c(NA,blockstats1k$time[-length(blockstats1k$time)])
blockstats1k$txpermin <- blockstats1k$transactions_count*60/blockstats1k$elapsed_time
blockstats1k$bytespermin <- blockstats1k$size*60/blockstats1k$elapsed_time

library(Hmisc)
describe(blockstats1k)
plot(blockstats1k$elapsed_time, blockstats1k$size, main="Block size v. time between blocks", xlab="Seconds since prior block discovery", ylab="bytes", pch=19)


#let's see the 100 largest blocks in our database
res <- dbSendQuery(con, "select height, size from blocks order by size desc limit 100")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
size <- dbFetch(res)
print(size)
dbClearResult(res)


#let's scatterplot all the blocks in our database by height and size
res <- dbSendQuery(con, "select height, size from blocks where branch=0")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
heightsize <- dbFetch(res)
dbClearResult(res)
plot(heightsize$height, heightsize$size, main="Bitcoin block sizes over time", 
   xlab="Block Height ", ylab="Size in Bytes ", pch=".")



#let's download all the data
res <- dbSendQuery(con, "select * from blocks where branch=0")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
blockchain <- dbFetch(res)
dbClearResult(res)
blockchain$work <- vapply(blockchain$work, paste, collapse = ", ", character(1L))
write.table(blockchain, file = "~/Desktop/blockchain.csv", sep = ",", col.names = NA,
            qmethod = "double")




dbDisconnect(con)
