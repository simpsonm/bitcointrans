## this is mostly a scratchpad right now
library(DBI)
library(rstan)

con <- dbConnect(RPostgres::Postgres(),dbname = 'toshi', 
                 host = 'toshi.cn6zzwcfsto5.us-east-1.rds.amazonaws.com',
                 port = 5432,
                 user = 'readonly',
                 password = 'password')

dbListTables(con)

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

save(timedata, file = "timedata.RData")
load("timedata.RData")

nobs <- 100
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

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=2000, chains=4, control = list(adapt_delta = 0.9, max_treedepth = 13)))

save(fit2, file="loggammafit.RData")

load("loggammafit.RData")

rstan::traceplot(fit2, pars=c('a', 'b', 'mu', 'sigma', 'lp__'), ncol=1)

print(fit2, pars=c('a', 'b', 'mu', 'sigma', 'lp__'), digits = 4)

fitex <- extract(fit2)



library(MCMCpack)

niter <- 4000
parout <- cbind(mu=fitex$mu, b=fitex$b, sigma=fitex$sigma)
parsum <- summary(mcmc(parout))

mus <- c(0, .01, .05, .1, .5)
K <- length(mus)
out <- matrix(0, ncol=3, nrow=K)


for(i in 1:K){
  ## use posterior means
  pars <- parsum[[1]][,1]
  muhat <- pars[1]
  bhat <- pars[2]
  sigmahat <- pars[3]
  xs <- rexp(niter, 1/10)
  taus <- cumsum(xs)
  mbs <- rgamma(niter, (muhat + mus[i])*bhat*xs, bhat)
  out[i,1] <- mean(mbs>1)

  ## use posterior medians
  pars <- parsum[[2]][,3]
  muhat <- pars[1]
  bhat <- pars[2]
  sigmahat <- pars[3]
  xs <- rexp(niter, 1/10)
  taus <- cumsum(xs)
  mbs <- rgamma(niter, (muhat + mus[i])*bhat*xs, bhat)
  out[i,2] <- mean(mbs>1)

  ## use full posterior
  xs <- rexp(niter, 1/10)
  taus <- cumsum(xs)
  mbs <- rgamma(niter, (fitex$mu + mus[i])*fitex$b*xs, fitex$b)
  out[i,3] <- mean(mbs>1)
}

colnames(out) <- c("Mean", "Median", "Posterior")
rownames(out) <- mus

library(xtable)

xtable(out)






omega <- 2*pi/(60*24)

gamlist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10, mb = mbs, sigpars = c(1, 1), sigepspars = c(1, 1), bpars=c(1,1), omega=omega, alphamns=rep(0,3), alphasds=rep(1,3), rhopars=c(1,1), phipar = 2*pi, lgam0pars = c(0, 1))
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'seasonallognormalgamma2.stan', data = gamlist, iter=1, chains=1)

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=100, chains=1, control = list(adapt_delta = 0.98, max_treedepth = 12)))

traceplot(fit2, pars=c('b', 'sigma', 'phi', 'sigeps', 'alpha0', 'alpha1', 'alpha2', 'rho'))

##, control = list(adapt_delta = 0.98, max_treedepth = 14)))

