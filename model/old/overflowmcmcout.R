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

lcgamma <- function(par, dat){
  k <- exp(par[1])
  alpha <- exp(par[2])
  beta <- exp(par[3])
  m <- dat[,1]
  d <- dat[,2]
  a <- m*k
  out <- - lbeta(alpha, a) - a*log(beta) + (a-1)*log(d) - (alpha + a)*log(1 + d/beta)
  return(-sum(out))
}


ahat <- optimize(lfun, c(0.00001, 1000), m = diff(datain$M), d = datain$D, maximum = TRUE)$maximum
bhat <- ahat/(sum(datain$D)/datain$M[nblocks + 1])

test <- optim(c(0,0,0), lcgamma, dat=cbind(diff(datain$M), datain$D))
testpars <- exp(test$par)
mmm <- tail(diff(datain$M), 1)
nnn <- 10000

outsim <- matrix(rgamma(nnn*mmm, rgamma(nnn*mmm, 1, 1.9), .07), nrow = mmm)
outsim <- matrix(rgamma(nnn*mmm, .04, .07), nrow = mmm)
outsim2 <- apply(outsim, 2, sum)

ssigma <- 1
vsigma <- 1
aphi <- 2
bphi <- 2
mugamma <- 3.5
sig2gamma <- 1
balpha <- 100 ##0000
aalpha <- 1*balpha
bbeta <- balpha
abeta <- (bhat/ahat*1)*bbeta
Dbars <- c(0.25, 0.35, 0.5, 0.75, 0.9, 0.95, 1)*1000

##Dbars <- c(1)*1000

priorin <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, mugamma = mugamma,
                sig2gamma = sig2gamma, aalpha = aalpha, balpha = balpha,
                abeta = abeta, bbeta = bbeta, aq = aalpha, bq = balpha, Dbars = Dbars)

initial <- list()
initial$tau <- datain$tobs[order(datain$tobs)] + cumsum(rexp(nblocks, 100000))
initial$sigma <- ssigma
initial$gam <- 0
initial$phi <- aphi/(aphi + bphi)
initial$N <- datain$M
initial$Dbar <- rep(max(priorin$Dbars), nblocks)
initial$alpha <- aalpha/balpha
initial$beta <- abeta/bbeta
initial$ps <- rep(ahat, tail(datain$M, 1) + 1)
initial$q <- bhat

niter <- 20000
rwc <- .1
tune <- TRUE

load("overout100k.RData")
draws <- test$draws
taucov <- cov(draws[,1:(nblocks) + 5 + nblocks])
tauchol <- t(chol(taucov))
rm(test)
rm(draws)

srwsds <- list(taulogrwsds = rep(0, nblocks), gamlogrwsd = -2.8, alphalogrwsd = -4.6,
               rwc = rwc, tune = tune, accprob = c(0.44, 0.44, 0.44))

frwsds <- list(taulogrwsd = -5.3, gamlogrwsd = -2.8, alphalogrwsd = -4.6,
               plogrwsds = rep(-4.6, tail(datain$M, 1) + 1), 
               tauchol = tauchol,
               rwc = rwc, tune = tune, accprob = c(0.44, 0.44, 0.44))

trendprior <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, mugamma = c(3,0),
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
                    tauchol = tauchol, rwc = rwc, tune = tune,
                    accprob = c(0.44, 0.44, 0.44, 0.44))

niter <- 10000

stautest <- staumcmcC(niter, datain, initial, priorin, srwsds)
overstautest <- overstaumcmcC(niter, datain, initial, priorin, srwsds)
ftautest <- ftaumcmcC(niter, datain, initial, priorin, frwsds)

overftautest <- overftaumcmcC(niter, datain, initial, priorin, frwsds)

ftautrendtest <- ftautrendmcmcC(niter, datain, trendinitial, trendprior, ftrendrwsds)

overftrendtest <- overftautrendmcmcC(niter, datain, trendinitial, trendprior, ftrendrwsds)

apply(stautest$tauaccs, 2, mean)
mean(stautest$gamaccs)
mean(stautest$alphaaccs)

apply(overstautest$tauaccs, 2, mean)
mean(overstautest$gamaccs)
mean(overstautest$alphaaccs)

mean(ftautest$tauaccs)
mean(ftautest$gamaccs)
mean(ftautest$alphaaccs)

mean(overftautest$tauaccs)
mean(overftautest$gamaccs)
mean(overftautest$alphaaccs)

mean(ftautrendtest$tauaccs)
mean(ftautrendtest$gamaccs)
mean(ftautrendtest$alphaaccs)

mean(overftrendtest$tauaccs)
mean(overftrendtest$gamaccs)
mean(overftrendtest$alphaaccs)

Ms <- datain$M[-1]
Ns <- overftautest$draws[,1:100 + 5]
NMs <- t(t(Ns) - Ms)

overprob <- apply(NMs[-c(1:5000),]>0, 2, mean)
idx <- which(overprob>0)

temp <- outer(Dbars, datain$D, '-')
temp[temp < 0] <- Inf

tempdat <- cbind(datain$D, apply(temp, 2, min))
tempdat <- cbind(tempdat, tempdat[,1] + tempdat[,2])
colnames(tempdat) <- c("size", "diff", "Dbar")
tempdat <- cbind(tempdat, overprob)



#### Dbar = block data cap
#### kb   = block data size
#### diff = difference between the two
#### based on my understanding of overflow, it only happens whenever there is a transaction
####  that is too large to fit in the block.
####  so any transaction that overflows has size >= diff when there is overflow.
####  that means there are transactions with, e.g., 45 kb of data.
####  would be great to get transaction size info to confirm this.




source("overflowmcmcCwrap.R")

rcpptest <- overtrendmcmcCwrap(niter, datain, trendinitial, trendprior, ftrendrwsds)




summary(mcmc(overftrendtest$draws[-c(1:5000),1:6]))[[1]]
summary(mcmc(rcpptest$draws[-c(1:5000),1:6]))[[1]]

summary(mcmc(overftrendtest$draws[-c(1:5000),1:6]))[[2]]
summary(mcmc(rcpptest$draws[-c(1:5000),1:6]))[[2]]


Ms <- datain$M[-1]
NsR <- overftrendtest$draws[,1:100 + 6]
NMsR <- t(t(NsR) - Ms)

overprobR <- apply(NMsR[-c(1:5000),]>0, 2, mean)
idxR <- which(overprobR>0)

NsC <- rcpptest$draws[,1:100 + 6]
NMsC <- t(t(NsC) - Ms)

overprobC <- apply(NMsC[-c(1:5000),]>0, 2, mean)
idxC <- which(overprobC>0)


tempR <- outer(Dbars, datain$D, '-')
tempR[tempR < 0] <- Inf

tempC <- outer(Dbars, datain$D, '-')
tempC[tempC < 0] <- Inf


tempdatR <- cbind(datain$D, apply(tempR, 2, min))
tempdatR <- cbind(tempdatR, tempdatR[,1] + tempdatR[,2])
colnames(tempdatR) <- c("size", "diff", "Dbar")
tempdatC <- cbind(datain$D, apply(tempC, 2, min))
tempdatC <- cbind(tempdatC, tempdatC[,1] + tempdatC[,2])
colnames(tempdatC) <- c("size", "diff", "Dbar")
tempdatR <- cbind(tempdatR, overprobR, (overprobR>0))
tempdatC <- cbind(tempdatC, overprobC, (overprobC>0))


tempdatR[tempdatR[,2]>5 & tempdatR[,2]<10,]


##### trying to see prior probability of overflow
npiter <- 1000
## sim pars
alpha <- rgamma(npiter, aalpha, balpha)
beta <- rgamma(npiter, abeta, bbeta)
gam1 <- rnorm(npiter, mugamma, sqrt(sig2gamma)) + 3.5
gam2 <- rnorm(npiter, mugamma, sqrt(sig2gamma))*0
omega <- rgamma(npiter, 1/2, 1/(2*vsigma*ssigma^2))
sigma2 <- 1/rgamma(npiter, vsigma/2, omega)
phi <- rbeta(npiter, aphi, bphi)
## now sim the rest
npblocks <- 100
Dbar <- matrix(sample(Dbars, npiter*npblocks, TRUE), ncol = npblocks)
delta <- matrix(rexp(npblocks*npiter, 1/10), ncol = npblocks)
tau <- t(apply(delta, 1, cumsum))
psi <- matrix(0, nrow = npiter, ncol = npblocks)
eta <- matrix(0, nrow = npiter, ncol = npblocks)
for(i in 1:npiter){
  psi[i,] <- psifunC(tau[i,], c(gam1[i], gam2[i]), c(0,tau[i,-length(tau[i,])]))
  eta[i,] <- rnbinom(npblocks, psi[i,], phi[i])
}
N <- t(apply(eta, 1, cumsum))

M <- matrix(0, nrow = npiter, ncol = npblocks)
D <- matrix(0, nrow = npiter, ncol = npblocks)
for(i in 1:npiter){
  print(i)
  d <- rgamma(N[i,ncol(N)], alpha[i], beta[i])
  idx <- 1
  for(k in 1:npblocks){
    Dsum <- 0
    while(Dsum + d[idx] <= Dbar[i,k] & idx < N[i,k]){
      Dsum <- Dsum + d[idx]
      idx <- idx + 1
    }
    D[i,k] <- Dsum
    M[i,k] <- idx
    idx <- idx + 1
  }
}


overprob <- apply((N - M)>0, 2, mean)
overidx <- which(N[,1] - M[,1] > 0)
cbind(D[overidx, 1], Dbar[overidx,1], Dbar[overidx,1] - D[overidx, 1])



source("overflowmcmc.R")

niter <- 10000

cgamtest <- cgammamcmcC(niter, datain,          initial, priorin,          frwsds)

cgamtest <- cgammamcmcC(niter, datain, cgamtest$initial, priorin, cgamtest$rwsds)

initialtemp <- overcgamtest$initial
rwsdstemp <- overcgamtest$rwsds
rm(overcgamtest)

rwsdstemp <- frwsds
frwsds$accprob <- c(0.234, 0.44, 0.44, 0.44)

overcgamtest <- overcgammamcmc(niter, datain, initialtemp, priorin, rwsdstemp)

kk <- 4 #nblocks + 3 - 4*0
par(mfcol=c(4,2))
coda::traceplot(mcmc(overftautest$draws[-c(1:5000),1:4+kk]))
coda::traceplot(mcmc(overcgamtest$draws[-c(1:5000),1:4+kk]))

par(mfcol=c(1,1))
coda::traceplot(mcmc(overcgamtest$draws[,7]))

kk <- Mtotal + 1 + 6
par(mfcol=c(3,2))
coda::traceplot(mcmc(overcgamtest$draws[,1:6 + kk]))

datain$D[
  which(apply(overcgamtest$draws[,kk + 1:nblocks], 2, var)>0)
]

kk <- nblocks - 1 - 6
par(mfcol=c(3,2))
coda::traceplot(mcmc(overftautest$draws[-c(1:3000),1:6+kk]))

datain$D[
  which(apply(overftautest$draws[,1:nblocks + 5], 2, var) > 0)
]



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


