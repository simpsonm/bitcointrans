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


# let's see if we can find the widely-used block size soft limits
res <- dbSendQuery(con, "select height, size from blocks where branch=0")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
heightsize <- dbFetch(res)
dbClearResult(res)
d<-density(heightsize$size)
plot(d)

# 5 widely-used block size soft limits: 250kb, 750kb, something in the 300kb range, something in the 900kb range, and 1mb (also occasionally 500kb)
# upon further research: 
# 250kb limit introduced in July 2012
# 350kb in version 0.8.6, released December 9, 2013 https://gist.github.com/gavinandresen/7670433#086-accept-into-block
# 750kb in version 0.9.0, released March 19, 2014 https://bitcoin.org/en/release/v0.9.0
# the value or values in the 900 range may be due almost exclusively to the Elegius mining pool. Still trying to pin down the exact value


# let's also see if we can find the stress tests programmatically:
# here's what we know
# stress test data
# 1. May 2015
# blocks became full starting at block number 358596
# remained full until block number 358609
# citation: http://www.ofnumbers.com/2015/05/31/a-few-results-from-the-first-intentional-intentional-stress-test-on-a-communal-blockchain/
# 2. June 22-23, 2015, 100 blocks
# citation: http://www.coindesk.com/bitcoin-network-survives-stress-test/
# 3. July 7-10, 2015
# https://medium.com/blockcypher-blog/a-bitcoin-spam-attack-post-mortem-s-la-ying-alive-654e914edcf4#.y2rcrl10n
# 4. September 2015
# citation: https://bitcoinmagazine.com/articles/coinwallet-crowdsources-transactions-major-stress-test-bitcoin-giveaway-1442011334

# install.packages("forecast")
library(forecast)
heightsize$mav <- ma(heightsize$size, 5000, centre=FALSE)
plot(heightsize$height, heightsize$size, main="Bitcoin block sizes over time", 
   xlab="Block Height ", ylab="Size in Bytes ", pch=".")
lines(heightsize$height, heightsize$mav, type='l', col="blue")
# produces a weird blue area toward the end of the series, am i doing something wrong?

# install.packages("zoo")
library(zoo)
heightsize$rollmed <- rollmedian(heightsize$size, 999, fill=NA, align="right")
plot(heightsize$height, heightsize$size, main="Bitcoin block sizes over time", 
   xlab="Block Height ", ylab="Size in Bytes ", pch=".")
lines(heightsize$height, heightsize$rollmedian, type='l', col="blue")



dbDisconnect(con)
