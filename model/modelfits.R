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

res <- dbSendQuery(con, "select branch, height, time, size, transactions_count from blocks where branch = 0 order by height desc limit 5011")
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

nobs <- 500
blockdatashort <- blockdata[(1000 - nobs + 1):1011,]
lbs <- rep(0,nobs)
for(i in 1:nobs){
  lbs[i] <- median(blockdatashort$min[i:(i+10)])
}
times <- blockdatashort$min[-c(1:11)]

sum(lbs > times) ## check that each lower bound is below what it's bounding, should be 0

times <- times - blockdatashort$min[11]
lbs <- lbs - blockdatashort$min[11]
mbs <- blockdatashort$size[-c(1:11)]/1000000  ## convert to megabytes

gamlist <- list(nobs = nobs, t = times, lb = lbs, lambda = 1/10, x = mbs, sigpars = c(0, 2.5), loggammapars = c(0, 10), betapars = c(1, 1))
## lambda = 1 / mean block discovery time = 1/10 minutes

fit <- stan(file = 'lognormalgamma.stan', data = gamlist, iter=1, chains=1)
## initializes the model and does some basic checks

system.time(fit2 <- stan(fit = fit, data = gamlist, iter=4000, chains=4, control = list(adapt_delta = 0.9, max_treedepth = 13)))

## not stored on github because it's large, csv file of just model pars below is.
save(fit2, file="loggammafit.RData")
load("loggammafit.RData") 

rstan::traceplot(fit2, pars=c('gamma', 'beta', 'sigma', 'lp__'), ncol=1)

print(fit2, pars=c('gamma', 'beta', 'sigma', 'lp__'), digits = 4)

## save just the model parameters so we don't have to keep fitting the model
fitex <- extract(fit2)
modelpar <- cbind(gamma=fitex$gamma, beta=fitex$beta, sigma=fitex$sigma, lp__=fitex$lp__)
write.csv(modelpar, file="postpar.csv", row.names = FALSE)

fitex <- read.csv("postpar.csv")

niter <- 8000
parout <- cbind(gamma=fitex$gamma, beta=fitex$beta, sigma=fitex$sigma)
parsum <- summary(mcmc(parout))
parout <- cbind(parsum[[1]][,c(1,2,4)], parsum[[2]])
rownames(parout) <- c("$\\gamma$", "$\\beta$", "$\\sigma$")
colnames(parout)[3] <- "SE(Mean)"

print(xtable(parout, digits = 4), sanitize.rownames.function = function(x)gsub("\\\\", "\\", x, fixed=TRUE))

gams <- seq(0.02, 0.09, 0.01)
betamed <- parout[2,6]
betamn <- parout[2,1]
gammn <- parout[1,1]
K <- length(gams)
out <- matrix(0, ncol=K, nrow=5)
nrep <- 100

for(i in 1:K){
  ## constant transaction rate
  out[1,i] <- exp(-1/10/gams[i])
  ## use posterior median * 3
  xs <- rexp(nrep*niter, 1/10)
  mbs <- rgamma(nrep*niter, gams[i]*betamed*xs*3, betamed*3)
  out[2,i] <- mean(mbs>1)
  ## use posterior median
  xs <- rexp(nrep*niter, 1/10)
  mbs <- rgamma(nrep*niter, gams[i]*betamed*xs, betamed)
  out[3,i] <- mean(mbs>1)
  ## use full posterior for beta
  xs <- rexp(nrep*niter, 1/10)
  mbs <- rgamma(nrep*niter, gams[i]*rep(fitex$beta,nrep)*xs, rep(fitex$beta,nrep))
  out[4,i] <- mean(mbs>1)
  ## use full posterior for both
  xs <- rexp(nrep*niter, 1/10)
  mbs <- rgamma(nrep*niter, (rep(fitex$gamma,nrep) + gams[i] - gammn)*rep(fitex$beta,nrep)*xs, rep(fitex$beta,nrep))
  out[5,i] <- mean(mbs>1)
}

colnames(out) <- gams
rownames(out) <- c("$\\beta = \\infty$", "$\\beta = 3\\times\\beta^{(0.5)}$", "$\\beta = \\beta^{(0.5)}$", "$\\beta \\sim p(\\beta|x,t)$", "$(\\beta,\\gamma) \\sim p(\\beta,\\gamma|x,t)$")

print(xtable(out, digits = 4), sanitize.rownames.function = function(x)gsub("\\\\", "\\", x, fixed=TRUE))

