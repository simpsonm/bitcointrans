library(DBI)
library(rstan)

con <- dbConnect(RPostgres::Postgres(),dbname = 'toshi', 
                 host = 'toshi.cn6zzwcfsto5.us-east-1.rds.amazonaws.com',
                 port = 5432,
                 user = 'readonly',
                 password = 'password')

res <- dbSendQuery(con, "select height, time, size, transactions_count from blocks order by height desc limit 1011")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}

timedata <- dbFetch(res)
dbClearResult(res)
timedata <- timedata[order(timedata$height),]
timedata$elapsed_time <- timedata$time - c(NA,timedata$time[-length(timedata$time)])
timedata$min <- timedata$time/60 ## convert to minutes


nobs <- 500 
timedatashort <- timedata[(1000 - nobs + 1):1011,]
lbs <- rep(0,nobs)
for(i in 1:nobs){
  lbs[i] <- median(timedatashort$min[i:(i+10)])
}
times <- timedatashort$min[-c(1:11)]

which(lbs > times) ## check that each lower bound is below what it's bounding

times <- times - timedatashort$min[11]
lbs <- lbs - timedatashort$min[11]

timelist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10)
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'time.stan', data = timelist, iter=1, chains=1)
system.time(fit2 <- stan(fit = fit, data = timelist, iter=1000, chains=2))


fit2out <- extract(fit2)
library(MCMCpack)

par(mfrow=c(4,4))
traceplot(mcmc(fit2out$sigma))
traceplot(mcmc(fit2out$tau[,1:15]))

summary(mcmc(fit2out$tau))

#### lognormal seems to work fine, need to fit bigger model, also investigate other error
#### distributions, e.g. gamma

