library(DBI)
library(rstan)
library(MCMCpack)
library(xtable)

con <- dbConnect(RPostgres::Postgres(),dbname = 'toshi', 
                 host = 'toshi.cn6zzwcfsto5.us-east-1.rds.amazonaws.com',
                 port = 5432,
                 user = 'readonly',
                 password = 'password')

dbListTables(con)

res <- dbSendQuery(con, "select height, time, size, transactions_count from blocks order by height desc limit 5011")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}

blockdata <- dbFetch(res)
dbClearResult(res)
blockdata <- blockdata[order(blockdata$height),]
blockdata$elapsed_time <- blockdata$time - c(NA,blockdata$time[-length(blockdata$time)])
blockdata$min <- blockdata$time/60 ## convert to minutes

save(blockdata, file = "blockdata.RData")

load("blockdata.RData")


############# Gamma fits

## fit on all 5000
mbs <- blockdata$size[-c(1:11)]/1000000  ## convert to megabytes
gamma5000list <- list(mb = mbs, logalphapars = c(0,10), betapars = c(1, 1), nobs = length(mbs))
g5000fit <- stan(file = 'gamma.stan', data = gamma5000list, chains = 4, iter = 4000)
g5000out <- extract(g5000fit)
g5000line <- c(summary(mcmc(matrix(g5000out$alpha, ncol=1)))[[1]][1], summary(mcmc(matrix(g5000out$alpha, ncol=1)))[[2]][c(1,5)])

## fit on blocks of 1000
g1000line <- NULL
for(i in 1:5){
  mbs <- blockdata$size[-c(1:11)]/1000000  ## convert to megabytes
  mbs <- mbs[1:1000 + 1000*(i-1)]
  gamma1000list <- list(mb = mbs, logalphapars = c(0,10), betapars = c(1, 1), nobs = length(mbs))
  g1000fit <- stan(file = 'gamma.stan', data = gamma1000list, chains = 4, iter = 4000)
  print(g1000fit)
  g1000out <- extract(g1000fit)
  g1000line <- rbind(g1000line, c(summary(mcmc(matrix(g1000out$alpha, ncol=1)))[[1]][1], summary(mcmc(matrix(g1000out$alpha, ncol=1)))[[2]][c(1,5)]))
}

## fit on blocks of 500
g500line <- NULL
for(i in 1:10){
  mbs <- blockdata$size[-c(1:11)]/1000000  ## convert to megabytes
  mbs <- mbs[1:500 + 500*(i-1)]
  gamma500list <- list(mb = mbs, logalphapars = c(0,10), betapars = c(1, 1), nobs = length(mbs))
  g500fit <- stan(file = 'gamma.stan', data = gamma500list, chains = 4, iter = 4000)
  print(g500fit)
  g500out <- extract(g500fit)
  g500line <- rbind(g500line, c(summary(mcmc(matrix(g500out$alpha, ncol=1)))[[1]][1], summary(mcmc(matrix(g500out$alpha, ncol=1)))[[2]][c(1,5)]))
}

gsummary <- rbind(g5000line, g1000line, g500line)
rownames(gsummary) <- c("n=5000", paste("n=1000, g=", 1:5, sep=""), paste("n=500, g=", 1:10, sep=""))
xtable(gsummary)

############## Full model fits

nobs <- 1000
blockdatashort <- blockdata[(1000 - nobs + 1):1011,]
lbs <- rep(0,nobs)
for(i in 1:nobs){
  lbs[i] <- median(blockdatashort$min[i:(i+10)])
}
times <- blockdatashort$min[-c(1:11)]
mbs <- blockdata$size[-c(1:11)]/1000000  ## convert to megabytes

sum(lbs > times) ## check that each lower bound is below what it's bounding, should be 0

times <- times - blockdatashort$min[11]
lbs <- lbs - blockdatashort$min[11]
mbs <- blockdatashort$size[-c(1:11)]/1000000  ## convert to megabytes

gamlist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10, x = mbs, sigpars = c(0, 2.5), loggammapars = c(0, 10), betapars = c(1, 1))
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'lognormalgamma.stan', data = gamlist, iter=1, chains=1)
## initializes the model and does some basic checks

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=4000, chains=4, control = list(adapt_delta = 0.9, max_treedepth = 13)))

save(fit2, file="loggammafit.RData")

load("loggammafit.RData")

rstan::traceplot(fit2, pars=c('a', 'b', 'mu', 'sigma', 'lp__'), ncol=1)

print(fit2, pars=c('a', 'b', 'mu', 'sigma', 'lp__'), digits = 4)

fitex <- extract(fit2)

niter <- 8000
parout <- cbind(mu=fitex$mu, b=fitex$b, sigma=fitex$sigma)
parsum <- summary(mcmc(parout))
parout <- cbind(parsum[[1]][,c(1,2,4)], parsum[[2]])
rownames(parout) <- c("$\\gamma$", "$\\beta$", "$\\sigma$")
colnames(parout)[3] <- "SE of Mean"

xtable(parout, digits = 4)


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

