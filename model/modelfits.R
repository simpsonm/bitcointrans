## this is mostly a scratchpad right now


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


nobs <- 1000
timedatashort <- timedata[(1000 - nobs + 1):1011,]
lbs <- rep(0,nobs)
for(i in 1:nobs){
  lbs[i] <- median(timedatashort$min[i:(i+10)])
}
times <- timedatashort$min[-c(1:11)]

which(lbs > times) ## check that each lower bound is below what it's bounding

times <- times - timedatashort$min[11]
lbs <- lbs - timedatashort$min[11]
mbs <- timedatashort$size[-c(1:11)]/1000000  ## convert to megabytes

gamlist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10, mb = mbs, sigpars = c(1, 1), apars=c(1, 1), bpars=c(1,1))
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'lognormalgamma.stan', data = gamlist, iter=1, chains=1)
## initializes the model and does some basic checks

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=200, chains=1, control = list(adapt_delta = 0.9, max_treedepth = 13)))


gamgamlist <- list(nobs = nobs, t = times/60/24/7/52, lb = lbs/60/24/7/52, lambda = 1/10/60/24/7/52, mb = mbs, betapars = c(1, 1), apars=c(1, 1), bpars=c(1,1))

initlist <- list(list(x = rep(1, length(times)), a = 1, b = 1, beta = 1))

fit <- stan(file = 'gammagamma.stan', data = gamgamlist, iter = 1, chains = 1)

system.time(fit2 <- stan(fit = fit, data = gamgamlist, iter=1, chains=1, control = list(adapt_delta = 0.9, max_treedepth = 13)))

rstan::traceplot(fit2, pars=c('a', 'b', 'mu', 'beta'), ncol=1)

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=4000, chains=4))




fitex <- extract(fit2)

