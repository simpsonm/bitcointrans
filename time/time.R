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

timelist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10)
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'time.stan', data = timelist, iter=1, chains=1)
## initializes the model and does some basic checks

system.time(fit2 <- stan(fit = fit, data = timelist, iter=2000, chains=4))
## fits the model using Hamiltonian Monto Carlo. Took about 45 minute son my laptop.
## note: easy to parallilize this across chains.



## traceplots to check for convergence, looks good
library(MCMCpack)
fit2out <- extract(fit2)
par(mfrow=c(4,4))
traceplot(mcmc(fit2out$sigma))
traceplot(mcmc(fit2out$tau[,1:15]))

## summarizes all estimates and gives some diagnostics
## neff = effective sample size. Closer to actual samples size is better
## Rhat is a measure of convergence, closer to 1 is better
print(fit2, pars = c('sigma', 'tau'))

## collect and report tau estimates for blocks near negative elapsed times
taumeans <- summary(mcmc(fit2out$tau))[[1]][,1]
tauquants <- summary(mcmc(fit2out$tau))[[2]][,c(1,5)]
idxs <- which(diff(times)<0)
idxs <- c(idxs, idxs-1, idxs+1)
idxs <- idxs[order(idxs)]
negativetime <- cbind(times[idxs], taumeans[idxs], tauquants[idxs,])
rownames(negativetime) <- idxs
colnames(negativetime) <- c('t', 'tau', '2.5%', '97.5%')

round(negativetime, 2)
