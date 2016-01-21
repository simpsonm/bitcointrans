source("overflowmcmc.R")
library(MCMCpack)

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
aphi <- 1/50
bphi <- 1/100
agamma <- 1
bgamma <- 1
balpha <- 1/10
aalpha <- ahat*balpha ## 1/10
bbeta <- 1/10
abeta <- bhat*bbeta
Dbars <- c(0.25, 0.35, 0.5, 0.75, 0.9, 0.95, 1)*1000
priorin <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, agamma = agamma,
                bgamma = bgamma,
                aalpha = aalpha, balpha = balpha, abeta = abeta, bbeta = bbeta, Dbars = Dbars)

initial <- list()
initial$tau <- datain$tobs[order(datain$tobs)] + cumsum(rexp(nblocks, 100000))
initial$sigma <- ssigma
initial$gam <- 1
initial$phi <- aphi/bphi
initial$N <- datain$M
initial$Dbar <- rep(max(priorin$Dbars), nblocks)
initial$alphad <- aalpha/balpha
initial$betad <- abeta/bbeta
initial$lambda <- rep(1, nblocks)

niter <- 10000
rwc <- .1
H <- 50
hightarget <- rep(0.45,3)
lowtarget <- rep(0.4,3)
tune <- TRUE

load("overout100k.RData")
draws <- test$draws
taucov <- cov(draws[,1:(nblocks) + 5 + nblocks])
tauchol <- t(chol(taucov))
rm(test)
rm(draws)

srwsds <- list(taulogrwsds = rep(0, nblocks), gamlogrwsd = 0, alphadlogrwsd = 0,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)
frwsds <- list(taulogrwsd = 0, gamlogrwsd = 0, alphadlogrwsd = 0, tauchol = tauchol,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)

stauout <- list()
ftauout <- list()
overstauout <- list()
overftauout <- list()

stest <- staumcmc(niter, datain, initial, priorin, srwsds)
print("stest finished")
ftest <- ftaumcmc(niter, datain, initial, priorin, frwsds)
print("ftest finished")

hightarget <- c(.27,rep(0.45,2))
lowtarget <- c(.22,rep(0.4,2))
srwsds <- list(taulogrwsds = rep(0, nblocks), gamlogrwsd = 0, alphadlogrwsd = 0,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)
frwsds <- list(taulogrwsd = 0, gamlogrwsd = 0, alphadlogrwsd = 0, tauchol = tauchol,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)

overstest <- overstaumcmc(niter, datain, initial, priorin, srwsds)
print("overstest finished")
overftest <- overftaumcmc(niter, datain, initial, priorin, frwsds)
print("overftest finished")

trendprior <- list(ssigma = ssigma, vsigma = vsigma, aphi = aphi, bphi = bphi, mugamma = c(0,0),
                sig2gamma = c(1,1),
                aalpha = aalpha, balpha = balpha, abeta = abeta, bbeta = bbeta, Dbars = Dbars)

trendinitial <- list()
trendinitial$tau <- datain$tobs[order(datain$tobs)] + cumsum(rexp(nblocks, 100000))
trendinitial$sigma <- ssigma
trendinitial$gam <- c(0,0)
trendinitial$phi <- aphi/bphi
trendinitial$N <- datain$M
trendinitial$Dbar <- rep(max(priorin$Dbars), nblocks)
trendinitial$alphad <- aalpha/balpha
trendinitial$betad <- abeta/bbeta
trendinitial$lambda <- rep(1, nblocks)

ftrendrwsds <- list(taulogrwsd = 0, gamlogrwsds = c(0, 0), alphadlogrwsd = 0, tauchol = tauchol,
               rwc = rwc, H = H, tune = tune, lowtarget = lowtarget, hightarget = hightarget)

ftrendtest <- ftautrendmcmcC(niter, datain, trendinitial, trendprior, ftrendrwsds)

sdraws <- stest$draws
fdraws <- ftest$draws
oversdraws <- overstest$draws
overfdraws <- overftest$draws

mean(1-stest$gamrejs)
mean(1-stest$alphadrejs)
summary(apply(1-stest$taurejs, 2, mean))

mean(1-overstest$gamrejs)
mean(1-overstest$alphadrejs)
summary(apply(1-overstest$taurejs, 2, mean))

mean(1-ftest$gamrejs)
mean(1-ftest$alphadrejs)
mean(1-ftest$taurejs)

mean(1-overftest$gamrejs)
mean(1-overftest$alphadrejs)
mean(1-overftest$taurejs)

exp(stest$rwsds$gamlogrwsd)
exp(stest$rwsds$alphadlogrwsd)
summary(exp(stest$rwsds$taulogrwsds))
which.max(exp(stest$rwsds$taulogrwsds))

exp(overstest$rwsds$gamlogrwsd)
exp(overstest$rwsds$alphadlogrwsd)
summary(exp(overstest$rwsds$taulogrwsds))
which.max(exp(overstest$rwsds$taulogrwsds))

exp(ftest$rwsds$gamlogrwsd)
exp(ftest$rwsds$alphadlogrwsd)
exp(ftest$rwsds$taulogrwsd)

exp(overftest$rwsds$gamlogrwsd)
exp(overftest$rwsds$alphadlogrwsd)
exp(overftest$rwsds$taulogrwsd)


par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 ]))
traceplot(mcmc(fdraws[, 1:5 ]))
traceplot(mcmc(oversdraws[, 1:5 ]))
traceplot(mcmc(overfdraws[, 1:5 ]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[-c(1:5000), 1:5 ]))
traceplot(mcmc(fdraws[-c(1:5000), 1:5 ]))
traceplot(mcmc(oversdraws[-c(1:5000), 1:5 ]))
traceplot(mcmc(overfdraws[-c(1:5000), 1:5 ]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 + 5 ]))
traceplot(mcmc(fdraws[, 1:5 + 5 ]))
traceplot(mcmc(oversdraws[, 1:5 + 5 ]))
traceplot(mcmc(overfdraws[, 1:5 + 5 ]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 + 5 + nblocks]))
traceplot(mcmc(fdraws[, 1:5 + 5 + nblocks]))
traceplot(mcmc(oversdraws[, 1:5 + 5 + nblocks]))
traceplot(mcmc(overfdraws[, 1:5 + 5 + nblocks]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 + 5 + nblocks + 15]))
traceplot(mcmc(fdraws[, 1:5 + 5 + nblocks + 15]))
traceplot(mcmc(oversdraws[, 1:5 + 5 + nblocks + 15]))
traceplot(mcmc(overfdraws[, 1:5 + 5 + nblocks + 15]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 + 5 + nblocks*2]))
traceplot(mcmc(fdraws[, 1:5 + 5 + nblocks*2]))
traceplot(mcmc(oversdraws[, 1:5 + 5 + nblocks*2]))
traceplot(mcmc(overfdraws[, 1:5 + 5 + nblocks*2]))

par(mfcol=c(5,4))
traceplot(mcmc(sdraws[, 1:5 + 5 + nblocks*3]))
traceplot(mcmc(fdraws[, 1:5 + 5 + nblocks*3]))
traceplot(mcmc(oversdraws[, 1:5 + 5 + nblocks*3]))
traceplot(mcmc(overfdraws[, 1:5 + 5 + nblocks*3]))


cbind(summary(mcmc(sdraws[-c(1:5000), 1:5]))[[1]][,1],
summary(mcmc(fdraws[-c(1:5000), 1:5]))[[1]][,1],
summary(mcmc(oversdraws[-c(1:5000), 1:5]))[[1]][,1],
summary(mcmc(overfdraws[-c(1:5000), 1:5]))[[1]][,1])



cbind(
  summary(mcmc(sdraws[-c(1:5000), 4] / sdraws[-c(1:5000),5]))[[1]][1],
  summary(mcmc(fdraws[-c(1:5000), 4] / fdraws[-c(1:5000),5]))[[1]][1],
  summary(mcmc(oversdraws[-c(1:5000), 4] / oversdraws[-c(1:5000),5]))[[1]][1],
  summary(mcmc(overfdraws[-c(1:5000), 4] / overfdraws[-c(1:5000),5]))[[1]][1]
)

cbind(
  summary(mcmc(sdraws[-c(1:5000), 4] / sdraws[-c(1:5000),5]^2))[[1]][1],
  summary(mcmc(fdraws[-c(1:5000), 4] / fdraws[-c(1:5000),5]^2))[[1]][1],
  summary(mcmc(oversdraws[-c(1:5000), 4] / oversdraws[-c(1:5000),5]^2))[[1]][1],
  summary(mcmc(overfdraws[-c(1:5000), 4] / overfdraws[-c(1:5000),5]^2))[[1]][1]
)

library(Hmisc)

sum(datain$D)/datain$M[nblocks + 1]

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

c(ahat, bhat)
ahat/bhat
ahat/(bhat^2)





par(mfrow=c(3,3))
traceplot(mcmc(draws)[,5 + nblocks + 1:9])

par(mfrow=c(3,3))
traceplot(mcmc(draws)[,5 + 2*nblocks - 1:9 + 1])


par(mfrow=c(3,3))
traceplot(mcmc(draws)[, 3*nblocks + 5 + 0:8 + 1])

par(mfrow=c(3,3))
traceplot(mcmc(draws)[, 4*nblocks + 5 - 0:8])


diff(summary(mcmc(draws[,5 + 1:nblocks]))[[1]][,1])


Nmns <- summary(mcmc(draws[,5 + 1:nblocks]))[[1]][,1]
names(Nmns) <- NULL
Nmns <- c(0, Nmns)

taumns <- summary(mcmc(draws[,5 + nblocks + 1:nblocks]))[[1]][,1]
names(taumns) <- NULL


tar <- function(x, a, b, L){
  dbeta(x, a, b)/pbeta(L, a, b, lower.tail=FALSE)
}

prop <- function(x, a, b, L){
  dbeta((x - L)/(1 - L), a, 1)/(1-L)
}

a <- .2
b <- 10
L <- .5
curve(tar(x, a, b, L), lwd=2, xlim=c(L,1))
curve(prop(x, a, b, L), lwd=2, col="blue", add=TRUE)
pbeta(L, a, b)

## idea: if beta < 1, check prej - use naive method if not too high, otherwise use beta trick
## then for alpha < 1 and not log concave

## still need to fix the alpha/beta/d steps somehow - alpha and beta are plummeting
## maybe try assuming that we just observe the Ns and see how it fits



out <- 0
for(i in 1:100){
  out <- out + (min(rgamma(4000, ahat)) == 0)
}
out/100

out <- 0
for(i in 1:100){
 out <- out + (min(exp(log(rgamma(4000, ahat + 1)) + (1/ahat)*log(runif(4000)))) == 0)
}
out/100


trad <- rgamma(4000, ahat)

test <- rgamma(4000000, ahat + 1)
min(test)

test2 <- runif(4000000)
min(test2)

nsam <- 100000
test3 <- log(runif(nsam))/ahat + log(rgamma(nsam, ahat + 1))
c(min(test3), max(test3))


Dk <- 932720
nsam <- 10000
lo <- log(rgamma(nsam, ahat + 1)) + log(runif(nsam))/ahat
lx <- lo + log(Dk) - log(sum(exp(lo)))
x <- exp(lx)
min(x)
