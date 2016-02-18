source("overflowmcmc.R")
source("overflowmcmcCwrap.R")

load("../charts/fullblocks.RData")
cleanblocks <- blocks[,c("height", "time", "transactions_count", "size")]
cleanblocks$elapsed.min<- c(NA, diff(cleanblocks$time))/60

nblocks <- 100

block1k <- cleanblocks[(nrow(cleanblocks) - nblocks - 10):nrow(cleanblocks),]
block1k$tobs <- cumsum(block1k$elapsed.min)
block1k$ms <- 0
for(i in 1:nblocks){
  block1k$ms[i + 11] <- median(block1k$tobs[i:(i + 10)])
}


datain <- list()
datain$tobs <- block1k$tobs[11 + 1:nblocks] - block1k$tobs[10]
datain$m <- block1k$ms[11 + 1:nblocks] - block1k$tobs[10]
datain$M <- c(0, cumsum(block1k$transactions_count[11 + 1:nblocks]))
datain$D <- block1k$size[11 + 1:nblocks]/1000 ## convert to kilobytes

lfun <- function(a, m, d){
  out <- rep(0, length(a))
  M <- sum(m)
  dbar <- sum(d)/M
  b <- a/dbar
  for(i in 1:length(a)){
    out[i] <- sum(dgamma(d, a[i]*m, b[i], log=TRUE))
  }
  return(out)
}

ahat <- optimize(lfun, c(0.00001, 1000), m = diff(datain$M), d = datain$D, maximum = TRUE)$maximum
bhat <- ahat/(sum(datain$D)/datain$M[nblocks + 1])

ssigma <- 1
vsigma <- 1
aphi <- 2
bphi <- 2
mugamma <- 0
sig2gamma <- 1
balpha <- 1000000
aalpha <- 1*balpha
bbeta <- balpha
abeta <- (bhat/ahat*1)*bbeta
Dbars <- c(0.25, 0.35, 0.5, 0.75, 0.9, 0.95, 1)*1000

Dbars <- c(1)*1000

priorin <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, mugamma = mugamma,
                sig2gamma = sig2gamma,
                aalpha = aalpha, balpha = balpha, abeta = abeta, bbeta = bbeta, Dbars = Dbars)

initial <- list()
initial$tau <- datain$tobs[order(datain$tobs)] + cumsum(rexp(nblocks, 100000))
initial$sigma <- ssigma
initial$gam <- 0
initial$phi <- aphi/(aphi + bphi)
initial$N <- datain$M
initial$Dbar <- rep(max(priorin$Dbars), nblocks)
initial$alpha <- aalpha/balpha
initial$beta <- abeta/bbeta

niter <- 20000
rwc <- .1
H <- 50
hightarget <- c(.27,rep(0.45,2))
lowtarget <- c(.22,rep(0.4,2))
tune <- TRUE

load("overout100k.RData")
draws <- test$draws
taucov <- cov(draws[,1:(nblocks) + 5 + nblocks])
tauchol <- t(chol(taucov))
rm(test)
rm(draws)

srwsds <- list(taulogrwsds = rep(0, nblocks), gamlogrwsd = -2.8, alphalogrwsd = -4.6,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)
frwsds <- list(taulogrwsd = -5.3, gamlogrwsd = -2.8, alphalogrwsd = -4.6, tauchol = tauchol,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)

trendprior <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, mugamma = c(0,0),
                sig2gamma = c(1,1),
                aalpha = aalpha, balpha = balpha, abeta = abeta, bbeta = bbeta, Dbars = Dbars)

trendinitial <- list()
trendinitial$tau <- datain$tobs[order(datain$tobs)] + cumsum(rexp(nblocks, 100000))
trendinitial$sigma <- ssigma
trendinitial$gam <- c(0,0)
trendinitial$phi <- aphi/(aphi + bphi)
trendinitial$N <- datain$M
trendinitial$Dbar <- rep(max(priorin$Dbars), nblocks)
trendinitial$alpha <- aalpha/balpha
trendinitial$beta <- abeta/bbeta

ftrendrwsds <- list(taulogrwsd = -5.3, gamlogrwsds = c(-2.8, -8.8), alphalogrwsd = -4.6,
                    tauchol = tauchol, rwc = rwc, H = H, tune = tune,
                    lowtarget = lowtarget, hightarget = hightarget)

niter <- 100

stautest <- staumcmcC(niter, datain, initial, priorin, srwsds)
overstautest <- overstaumcmcC(niter, datain, initial, priorin, srwsds)
ftautest <- ftaumcmcC(niter, datain, initial, priorin, frwsds)
overftautest <- overftaumcmcC(niter, datain, initial, priorin, frwsds)
ftautrendtest <- ftautrendmcmcC(niter, datain, trendinitial, trendprior, ftrendrwsds)
overftrendtest <- overftautrendmcmcC(niter, datain, trendinitial, trendprior, ftrendrwsds)
rcpptest <- overtrendmcmcCwrap(niter, datain, trendinitial, trendprior, ftrendrwsds)


k <- nblocks + 5
par(mfrow=c(3,3))
coda::traceplot(mcmc(stautest$draws)[,k])
coda::traceplot(mcmc(overstautest$draws)[,k])
coda::traceplot(mcmc(ftautest$draws)[,k])
coda::traceplot(mcmc(overftautest$draws)[,k])
kk <- ifelse(k < 4, k, k + 1)
coda::traceplot(mcmc(ftautrendtest$draws)[,kk])
coda::traceplot(mcmc(overftrendtest$draws)[,kk])
coda::traceplot(mcmc(rcpptest$draws)[,kk])


