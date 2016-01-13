library(MCMCpack)
library(truncnorm) ## rtruncnorm is faster than rtnorm in the msm package

## steps left: tau, beta, alpha_d
## all will be ARS / ARMS / Metropolis steps
## then put it all together and test

aomega <- (vsigma + 1)/2
bomega <- 1/(vsigma + s^2)

tfull <- list()

blocks <- matrix(1:nblocks, ncol = 1)

listdiffsq <- function(k, tlist, tauvec){
  out <- sum((tlist[[k]] - tauvec[k])^2)
  return(out)
}

## DA step: \tilde{n}_k and \tilde{t}_k
## Easily parallelizable
bsigma <- 0 ## running tally, used in sigma step
for(k in 1:nblock){
  rho <- pnorm((m[k] - tau[k])/sigma) ## rho = P(t_k < m_k) = P(fail)
  nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
  if(nrej[k] > 0){
    tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
    tfull[[k]] <- c(tmiss, tobs[k])
  } else {
    tfull[[k]] <- tobs[k]
  }
  bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
}

## sigma / sigma^2 step
omega <- rgamma(aomega, bomega + 1/sigma)
asigma <- (nblocks + sum(nrej) + vsigma)/2
sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
sigma <- sqrt(sigma^2)

## tau step -> complicated, probably not parallelizable

## phi step
phi <- rgamma(1, aphi + sum(psi), bphi + sum(lambda))

## beta step -> complicated

## d step
nds <- diff(M) ## M[1] = 0 and N[1] = 0
d <- list()
for(k in 1:nblocks){
  if(N[k] == M[k]){
    d[[k]] <- rdirichlet(1, rep(alphad, nds[k]))*D[k]
  } else {
    u <- runif(1)
    pupper <- pbeta(d[k], alphad, betad)
    d1 <- qbeta(pupper + u * (1 - pupper), alphad, betad)
    dother <- rdirichlet(1, rep(alphad, nds[k] - 1))
    d[[k]] <- c(d1, (1 - d1)*dother)*D[k]
  }
}
ds <- unlist(d) ## note to self: works exact as expected

## N step
for(k in 1:nblocks){
  if(sum(ds[(M[k] + 1):M[k+1]]) <= Dcap[k]){
    N[k+1] = M[k+1]
  } else {
    lower <- M[k+1] - N[k] 
    plower <- ppois(lower - 1, lambda[k])
    u <- runif(1)
    N[k+1] <- N[k] + qpois(u * (1 - plower) + plower, lambda[k])
  }
}


## Dbar step

## alphad step -> complicated

## betad step
betad <- rgamma(alphad*N[nblocks + 1] + abeta, bbeta + sum(d))

