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


dbDisconnect(con)