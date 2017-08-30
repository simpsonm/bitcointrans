library(MCMCpack)
library(truncnorm) ## rtruncnorm() is faster than rtnorm() in the msm package
library(ars)
library(compiler)
library(aster)

## overflow compound gamma testing
overcgammamcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  print("initial")
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0, tau))
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  ps <- initial$ps
  q <- initial$q
  psi <- psifunC(tau, c(gam,0), tau0)
  print("rwsds")
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  plogrwsds <- rwsds$plogrwsds  
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  prwsds <- exp(plogrwsds)  
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  print("probably expensive")
  paccs <- matrix(0, ncol = length(initial$ps), nrow = niter)
  print("expensive done")
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  tauchol <- rwsds$tauchol
  print("priors")
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aq <- prior$aq
  bq <- prior$bq
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  print("data")
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  rateq <- bq + sum(D)
  print("store stuff")
  ## Store draws: tau, sigma, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = 1 + nblocks + 1 + 1 + nblocks + nblocks + 2 + 1 + Mtotal + 1)
  print("expensive done")
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  pnames <- paste("p", 0:Mtotal + 1, sep="")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta", "q", pnames,
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  print("initializing DA stuff")
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  print("loop")
  for(iter in 1:niter){
    ## print(iter)
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    psiprop <- psifunC(tauprop, c(gam,0), tau0prop)
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    probtau <- exp(min(lanum - ladenom, 0))
    u <- runif(1)
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }
    ## print("tau step finished")
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }
    ## print("gamma step finished")
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(nds[k] > 1){
        pk <- ps[(M[k] + 1):M[k+1]]
        if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
          lo <- rloggammaC(nds[k], pk, 1)
          lo <- lo - max(lo)
          ld <- log(D[k]) + lo - log(sum(exp(lo)))
          ldlist[[k]] <- ld
        } else {
          lowlim <- (Dbar[k-1] - D[k-1])/D[k]
          d1 <- rbetarejC(1, pk[1], sum(pk[-1]), lowlim)
          loother <- rloggammaC(nds[k] - 1, pk[-1], 1)
          ldother <- loother - log(sum(exp(loother)))
          ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
          ldlist[[k]] <- ld
        }
      } else {
        ldlist[[k]] <- log(D[k])
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    Nover <- Ntotal - Mtotal
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, ps[Mtotal + 1], q)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, ps[Mtotal + 1], q, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), ps[Mtotal + 1], q))
      }
    } else {
      ldover <- rloggammaC(1, ps[Mtotal + 1], q)
    }
    ldlist[[nblocks + 1]] <- ldover
    ldvec <- c(ldvec, ldover)
    ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
      } else {
        lower <- max(M[k+1], N[k]) - N[k]
        psik <- psifunC(tau[k], c(gam, 0), tau0[k])
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          psik1 <- psifunC(tau[k+1], c(gam, 0), tau0[k+1])
          etaset <- lower:upper
          etalogprobs <- lgamma(etaset + psik) - lfactorial(etaset) +
            lgamma(upper - etaset + psik1) - lfactorial(upper - etaset)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          if(lower > 0){
            etak <- rktnb(1, psik, lower - 1, psik*(1-phi)/phi)
          } else {
            etak <- rnbinom(1, psik, phi)
          }
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    eta <- diff(N)
    ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
        Dbar[k] <- Dbars[min(idxk)]
      } else {
        idxk <- which(Dbars >= D[k])
        Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
      }
    }
    ## print("Dbar step finished")
    ## p step
    pprops <- exp(log(ps) + rnorm(Mtotal + 1, 0, prwsds))
    u <- runif(Mtotal + 1)
    pcons <- log(q) - beta + ldvec
    laps <- (pprops - ps)*pcons + alpha*log(pprops/ps) - lgamma(pprops) + lgamma(ps)
    probp <- exp(apply(laps, 2, function(x) max(x,0)))
    idxacc <- u < probp
    ps[idxacc] <- pprops[idxacc]
    paccs[iter, idxacc] <- 1
    ## print("p step finished")
    ## q step
    psum <- sum(ps)
    q <- rgamma(1, aq + psum, rateq + exp(ldover))
    ## beta step
    beta <- rgamma(1, abeta + alpha*(Mtotal + 1), bbeta + psum)
    ## print("q & beta steps finished")
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    con1 <- (Mtotal + 1)*log(beta) + sum(log(ps)) - balpha
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      (Mtotal + 1)*(lgamma(alphaprop) - lgamma(alpha))
    probalpha <- exp(min(laalpha, 0))
    u <- runif(1)
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }
    ## print("alpha step finished")
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
      plogrwsds <- plogrwsds + rwc*(probp - accprob[4])
      prwsds <- exp(plogrwsds)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, q, ps, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta, q = q, ps = ps,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd = taulogrwsd
  rwsds$gamlogrwsd = gamlogrwsd
  rwsds$alphalogrwsd = alphalogrwsd
  rwsds$plogrwsds = plogrwsds
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs, paccs = paccs)
  return(out)
}

## compound gamma testing
cgammamcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0,tau))
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  ps <- initial$ps
  q <- initial$q
  psi <- psifunC(tau, c(gam, 0), tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  plogrwsd <- rwsds$plogrwsd
  prwsd <- exp(plogrwsd)
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  print("probably expensive")
  paccs <- matrix(0, ncol = 1, nrow = niter)
  print("expensive done")
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aq <- prior$aq
  bq <- prior$bq
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  rateq <- bq + sum(D)
  ## Store draws: tau, sigma, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = 1 + nblocks + 1 + 1 + nblocks + nblocks + 2 + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta", "p", "q",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## print(iter)
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    psiprop <- psifunC(tauprop, c(gam, 0), tau0prop)
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    u <- runif(1)
    probtau <- exp(min(lanum - ladenom, 0))
    u <- runif(1)
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }    
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(nds[k] > 1){
        if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
          lo <- rloggammaC(nds[k], p, 1)
          lo <- lo - max(lo)
          ld <- log(D[k]) + lo - log(sum(exp(lo)))
          ldlist[[k]] <- ld
        } else {
          lowlim <- (Dbar[k-1] - D[k-1])/D[k]
          d1 <- rbetarejC(1, p, p*(nds[k] - 1), lowlim)
          loother <- rloggammaC(nds[k] - 1, p, 1)
          ldother <- loother - log(sum(exp(loother)))
          ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
          ldlist[[k]] <- ld
        }
      } else {
        ldlist[[k]] <- log(D[k])
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    ## p step
    pprop <- exp(log(p) + rnorm(1, 0, prwsd))
    u <- runif(1)
    con1 <- -beta + Ntotal*log(q) + Ntotal*mean(ldvec)
    lap <- (pprop - p)*con1 + (alpha - 1)*log(pprop/p) -
      Ntotal*(lgamma(pprop) - lgamma(p))
    probp <- exp(min(lap, 0))
    if(u < probp){
      p <- pprop
      paccs[iter,1] <- 1
    }
    ## q step
    q <- rgamma(1, aq + p*Ntotal, rateq)
    ## beta step
    beta <- rgamma(1, abeta + alpha, bbeta + p)
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    con1 <- -balpha + log(beta) + log(p)
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      (lgamma(alphaprop) - lgamma(alpha))
    u <- runif(1)
    probalpha <- exp(min(laalpha, 0))
    u <- runif(1)
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }
    ## Tuning update step
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
      plogrwsd <- plogrwsd + rwc*(probp - accprob[4])
      prwsd <- exp(plogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, p, q, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta, p = p, q = q,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd = taulogrwsd
  rwsds$gamlogrwsd = gamlogrwsd
  rwsds$alphalogrwsd = alphalogrwsd
  rwsds$plogrwsd = plogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs, paccs = paccs)
  return(out)
}


## no trend, no overflow, full tau step
ftaumcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0,tau))
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, c(gam, 0), tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  tauchol <- rwsds$tauchol
  accprob <- rwsds$accprob
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ratebeta <- bbeta + sum(D)
  ## Store draws: tau, sigma, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = 1 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## print(iter)
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    psiprop <- psifunC(tauprop, c(gam, 0), tau0prop)
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    probtau <- exp(min(lanum - ladenom, 0))
    u <- runif(1)
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }    
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- -balpha + Ntotal*log(beta) + sum(eta*log(D))
    laalpha <- (alphaprop - alpha)*con1 + (aalpha)*log(alphaprop/alpha) -
      (sum(lgamma(eta*alphaprop)) - sum(lgamma(eta*alpha)))
    probalpha <- exp(min(laalpha, 0))
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }    
    ## beta step
    beta <- rgamma(1, abeta + alpha*Ntotal, ratebeta)
    ## Tuning update step
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd <- taulogrwsd
  rwsds$gamlogrwsd <- gamlogrwsd
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

## no trend, overflow, full f step
overftaumcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0, tau))
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, c(gam,0), tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ratebeta <- bbeta + sum(D)
  ## Store draws: tau, sigma, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## print(iter)
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)/2
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    psiprop <- psifunC(tauprop, c(gam,0), tau0prop)
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    probtau <- exp(min(lanum - ladenom, 0))
    u <- runif(1)
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }    
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }    
    ## print("gamma step finished")
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(nds[k] > 1){
        if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
          lo <- rloggammaC(nds[k], alpha, 1)
          lo <- lo - max(lo)
          ld <- log(D[k]) + lo - log(sum(exp(lo)))
          ldlist[[k]] <- ld
        } else {
          lowlim <- (Dbar[k-1] - D[k-1])/D[k]
          d1 <- rbetarejC(1, alpha, alpha*(nds[k] - 1), lowlim)
          loother <- rloggammaC(nds[k] - 1, alpha, 1)
          ldother <- loother - log(sum(exp(loother)))
          ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
          ldlist[[k]] <- ld
        }
      } else {
        ldlist[[k]] <- log(D[k])
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    Nover <- Ntotal - Mtotal
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, alpha, beta)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, alpha, beta, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), alpha, beta))
      }
    } else {
      ldover <- rloggammaC(1, alpha, beta)
    }
    ldlist[[nblocks + 1]] <- ldover
    ldvec <- c(ldvec, ldover)
    ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
      } else {
        lower <- max(M[k+1], N[k]) - N[k]
        psik <- psifunC(tau[k], c(gam, 0), tau0[k])
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          psik1 <- psifunC(tau[k+1], c(gam, 0), tau0[k+1])
          etaset <- lower:upper
          etalogprobs <- lgamma(etaset + psik) - lfactorial(etaset) +
            lgamma(upper - etaset + psik1) - lfactorial(upper - etaset)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          if(lower > 0){
            etak <- rktnb(1, psik, lower - 1, psik*(1-phi)/phi)
          } else {
            etak <- rnbinom(1, psik, phi)
          }
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    eta <- diff(N)
    ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
        Dbar[k] <- Dbars[min(idxk)]
      } else {
        idxk <- which(Dbars >= D[k])
        Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
      }
    }
    ## print("Dbar step finished")
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- Ntotal*log(beta) + Ntotal*mean(ldvec) - balpha
    ## print(c(mean(ldvec), mean(exp(ldvec))))
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      Ntotal*(lgamma(alphaprop) - lgamma(alpha))
    probalpha <- exp(min(laalpha, 0))
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }        
    ## beta step
    beta <- rgamma(1, abeta + alpha*Ntotal, ratebeta + exp(ldover))
    ## Tuning update step
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd <- taulogrwsd
  rwsds$gamlogrwsd <- gamlogrwsd
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

## no trend, overflow, single tau steps
overstaumcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0, tau))
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, c(gam, 0), tau0)
  ## RW sds and bookkeeping
  taulogrwsds <- rwsds$taulogrwsds
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsds <- exp(taulogrwsds)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = nblocks, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    for(k in 1:nblocks){
      deltakold <- delta[k]
      zetakold <- log(deltakold)
      tausd <- taurwsds[k]
      zetakprop <- rnorm(1, zetakold, tausd)
      deltakprop <- exp(zetakprop)
      deltaprop <- delta
      deltaprop[k] <- deltakprop
      tauprop <- cumsum(deltaprop)
      tau0prop <- c(0, tauprop[-nblocks])
      psiprop <- psifunC(tauprop, c(gam, 0), tau0prop)
      oldsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tau))
      newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
      lanum <-   -newsqdiff/(2*sigma2) + psiprop[k]*log(phi) - deltakprop/10 +
        lgamma(psiprop[k] + eta[k]) - lgamma(psiprop[k]) + zetakprop
      ladenom <- -oldsqdiff/(2*sigma2) + psi[k]*log(phi)     - deltakold/10 +
        lgamma(psi[k] + eta[k])     - lgamma(psi[k])     + zetakold
      u <- runif(1)
      probtau <- exp(min(lanum - ladenom, 0))
      if(u < probtau){
        delta[k] <- deltakprop
        tau <- tauprop
        psi <- psiprop
        tau0 <- tau0prop
        tauaccs[iter, k] <- 1
      }
      if(tune){
        taulogrwsds[k] <- taulogrwsds[k] + rwc*(probtau - accprob[1])
        taurwsds[k] <- exp(taulogrwsds[k])
      }
    }
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposal
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }        
    ## print("gamma step finished")
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(nds[k] > 1){
        if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
          lo <- rloggammaC(nds[k], alpha, 1)
          lo <- lo - max(lo)
          ld <- log(D[k]) + lo - log(sum(exp(lo)))
          ldlist[[k]] <- ld
        } else {
          lowlim <- (Dbar[k-1] - D[k-1])/D[k]
          d1 <- rbetarejC(1, alpha, alpha*(nds[k] - 1), lowlim)
          loother <- rloggammaC(nds[k] - 1, alpha, 1)
          ldother <- loother - log(sum(exp(loother)))
          ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
          ldlist[[k]] <- ld
        }
      } else {
        ldlist[[k]] <- log(D[k])
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    Nover <- Ntotal - Mtotal
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, alpha, beta)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, alpha, beta, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), alpha, beta))
      }
      if(Nover > 1){
        ldover <- c(ldover, rloggammaC(Nover-1, alpha, beta))
      }
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    } else {
      ldover <- rloggammaC(1, alpha, beta)
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    }
    ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
      } else {
        lower <- max(M[k+1], N[k]) - N[k]
        psik <- psifunC(tau[k], c(gam,0), tau0[k])
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          psik1 <- psifunC(tau[k+1], c(gam,0), tau0[k+1])
          etaset <- lower:upper
          etalogprobs <- lgamma(etaset + psik) - lfactorial(etaset) +
            lgamma(upper - etaset + psik1) - lfactorial(upper - etaset)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          if(lower > 0){
            etak <- rktnb(1, psik, lower - 1, psik*(1-phi)/phi)
          } else {
            etak <- rnbinom(1, psik, phi)
          }
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    eta <- diff(N)
    ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
        Dbar[k] <- Dbars[min(idxk)]
      } else {
        idxk <- which(Dbars >= D[k])
        Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
      }
    }
    ## alpha step - Metrop
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- Ntotal*mean(ldvec) - Ntotal*log(bbeta + Ntotal*mean(exp(ldvec))) - balpha
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      Ntotal*(lgamma(alphaprop) - lgamma(alpha)) +
      lgamma(alphaprop*Ntotal + abeta - 1) - lgamma(alpha*Ntotal + abeta - 1)
    probalpha <- exp(min(laalpha, 0))
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }
    ## beta step
    ratebeta <- bbeta + sum(D)
    if(Nover > 0){
      ratebeta <- ratebeta + sum(exp(ldover))
    }
    beta <- rgamma(1, abeta + alpha*Ntotal, ratebeta)
    ## Tuning update step
    ## print("beta step finished")
    if(tune){
      ## tau tuned in its loop
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }    
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsds <- taulogrwsds
  rwsds$gamlogrwsd <- gamlogrwsd
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

## no trend, no overflow, single tau steps
staumcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  nblocks <- length(tau)
  tau0 <- c(0, tau[-nblocks])
  delta <- diff(c(0,tau))
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, c(gam, 0), tau0)
  ## RW sds and bookkeeping
  taulogrwsds <- rwsds$taulogrwsds
  gamlogrwsd <- rwsds$gamlogrwsd
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsds <- exp(taulogrwsds)
  gamrwsd <- exp(gamlogrwsd)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = nblocks, nrow = niter)
  gamaccs <- matrix(0, ncol = 1, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ratebeta <- bbeta + sum(D)
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
      nrej[k] <- rnbinom(1, 1, prob = 1 - rho)
      if(nrej[k] > 0){
        tmiss <- rtruncnorm(nrej[k], mean = tau[k], sd = sigma, a = -Inf, b = m[k])
        tfull[[k]] <- c(tmiss, tobs[k])
      } else {
        tfull[[k]] <- tobs[k]
      }
      bsigma <- bsigma + sum((tfull[[k]] - tau[k])^2)
    }
    ## print("DA Step Finished")
    ## sigma / sigma^2 step
    omega <- rgamma(aomega, bomega + 1/sigma)
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    for(k in 1:nblocks){
      deltakold <- delta[k]
      zetakold <- log(deltakold)
      tausd <- taurwsds[k]
      zetakprop <- rnorm(1, zetakold, tausd)
      deltakprop <- exp(zetakprop)
      deltaprop <- delta
      deltaprop[k] <- deltakprop
      tauprop <- cumsum(deltaprop)
      tau0prop <- c(0, tauprop[-nblocks])
      psiprop <- psifunC(tauprop, c(gam, 0), tau0prop)
      oldsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tau))
      newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
      lanum <-   -newsqdiff/(2*sigma2) + psiprop[k]*log(phi) - deltakprop/10 +
        lgamma(psiprop[k] + eta[k]) - lgamma(psiprop[k]) + zetakprop
      ladenom <- -oldsqdiff/(2*sigma2) + psi[k]*log(phi)     - deltakold/10 +
        lgamma(psi[k] + eta[k])     - lgamma(psi[k])     + zetakold
      u <- runif(1)
      probtau <- exp(min(lanum - ladenom, 0))
      if(u < probtau){
        delta[k] <- deltakprop
        tau <- tauprop
        psi <- psiprop
        tau0 <- tau0prop
        tauaccs[iter, k] <- 1
      }
      if(tune){
        taulogrwsds[k] <- taulogrwsds[k] + rwc*(probtau - accprob[1])
        taurwsds[k] <- exp(taulogrwsds[k])
      }      
    }
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    u <- runif(1)
    psiprop <- psifunC(tau, c(gamprop, 0), tau0)
    lanum <- sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) +
      sum(psiprop)*log(phi) - (gamprop - mugam)^2/sig2gam/2
    ladenom <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam - mugam)^2/sig2gam/2
    probgam <- exp(min(lanum - ladenom, 0))
    if(u < probgam){
      gam <- gamprop
      psi <- psiprop
      gamaccs[iter,1] <- 1
    }
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- -balpha + Ntotal*log(beta) + sum(eta*log(D))
    laalpha <- (alphaprop - alpha)*con1 + (aalpha)*log(alphaprop/alpha) -
      (sum(lgamma(eta*alphaprop)) - sum(lgamma(eta*alpha)))
    probalpha <- exp(min(laalpha, 0))
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }    
    ## beta step
    beta <- rgamma(1, abeta + alpha*Ntotal, ratebeta)
    ## Tuning update step
    if(tune){
      ## tau tuned in its loop
      gamlogrwsd <- gamlogrwsd + rwc*(probgam - accprob[2])
      gamrwsd <- exp(gamlogrwsd)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }    
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsds <- taulogrwsds
  rwsds$gamlogrwsd <- gamlogrwsd
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

## trend, no overflow, full tau step
ftautrendmcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  tau0 <- tau0[-(nblocks+1)]
  nblocks <- length(tau)
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, gam, tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsds <- rwsds$gamlogrwsds
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsds <- exp(gamlogrwsds)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 2, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = 2 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", paste("gam", 1:2, sep=""), "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
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
    asigma <- (nblocks + sum(nrej) + vsigma)
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    psiprop <- psifunC(tauprop, gam, tau0prop)
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    u <- runif(1)
    probtau <- exp(min(lanum - ladenom, 0))
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }    
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(2, gam, gamrwsds)
    u <- runif(2)
    psi1prop <- psifunC(tau, c(gamprop[1], gam[2]), tau0)
    lanum1 <- sum(lgamma(psi1prop + eta)) - sum(lgamma(psi1prop)) +
      sum(psi1prop)*log(phi) - (gamprop[1] - mugam[1])^2/sig2gam[1]/2
    ladenom1 <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam[1] - mugam[1])^2/sig2gam[1]/2
    probgam1 <- exp(min(lanum1 - ladenom1, 0))
    if(u[1] < probgam1){
      gam[1] <- gamprop[1]
      psi <- psi1prop
      gamaccs[iter,1] <- 1
    }        
    psi2prop <- psifunC(tau, c(gam[1], gamprop[2]), tau0)
    lanum2 <- sum(lgamma(psi2prop + eta)) - sum(lgamma(psi2prop)) +
      sum(psi2prop)*log(phi) - (gamprop[2] - mugam[2])^2/sig2gam[2]/2
    ladenom2 <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam[2] - mugam[2])^2/sig2gam[2]/2
    probgam2 <- exp(min(lanum2 - ladenom2, 0))
    if(u[2] < probgam2){
      gam[2] <- gamprop[2]
      psi <- psi2prop
      gamaccs[iter,2] <- 1
    }    
    lagam2 <- lanum2 - ladenom2
    ## alpha step - Metrop
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- -balpha + Ntotal*log(beta) + sum(eta*log(D))
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      sum(lgamma(eta*alphaprop) - lgamma(eta*alpha))
    probalpha <- exp(min(laalpha, 0))
    u <- runif(1)
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }
    ## print("alpha step finished")
    ## beta step
    beta <- rgamma(1, abeta + alpha*Ntotal, bbeta + sum(D))
    ## Tuning update step
    ## print("beta step finished")
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsds <- gamlogrwsds + rwc*(c(probgam1, probgam2) - accprob[2])
      gamrwsds <- exp(gamlogrwsds)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd <- taulogrwsd
  rwsds$gamlogrwsds <- gamlogrwsds
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

## trend, overflow, full tau step
overftautrendmcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  tau0 <- tau0[-(nblocks+1)]
  nblocks <- length(tau)
  sigma <- initial$sigma
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  eta <- diff(N)
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alpha <- initial$alpha
  beta <- initial$beta
  psi <- psifunC(tau, gam, tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsds <- rwsds$gamlogrwsds
  alphalogrwsd <- rwsds$alphalogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsds <- exp(gamlogrwsds)
  alpharwsd <- exp(alphalogrwsd)
  tauaccs <- matrix(0, ncol = 1, nrow = niter)
  gamaccs <- matrix(0, ncol = 2, nrow = niter)
  alphaaccs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  tune <- rwsds$tune
  accprob <- rwsds$accprob
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  mugam <- prior$mugamma
  sig2gam <- prior$sig2gamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(2*vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alpha, beta
  draws <- matrix(0, nrow = niter, ncol = 2 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", paste("gam", 1:2, sep=""), "alpha", "beta",
                       Nnames, taunames, Dbarnames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    ## print(iter)
    ## DA step: \tilde{n}_k and \tilde{t}_k
    ## Easily parallelizable
    bsigma <- 0 ## running tally, used in sigma step
    for(k in 1:nblocks){
      rho <- pnorm(m[k], tau[k], sigma) ## rho = P(t_k < m_k) = P(fail)
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
    asigma <- nblocks + sum(nrej) + vsigma
    sigma2 <- 1/rgamma(1, asigma/2, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    psiprop <- psifunC(tauprop, gam, tau0prop)
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    lanum <- -newsqdiff/(2*sigma2) + sum(psiprop)*log(phi) - tauprop[nblocks]/10 +
      sum(lgamma(psiprop + eta)) - sum(lgamma(psiprop)) + sum(zetaprop)
    ladenom <- -bsigma/(2*sigma2)+ sum(psi)*log(phi) - tau[nblocks]/10 +
      sum(lgamma(psi + eta)) - sum(lgamma(psi)) + sum(zetaold)
    u <- runif(1)
    probtau <- exp(min(lanum - ladenom, 0))
    if(u < probtau){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
      tauaccs[iter] <- 1
    }
    ## phi step
    phi <- rbeta(1, aphi + sum(psi), bphi + Ntotal)
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(2, gam, gamrwsds)
    u <- runif(2)
    psi1prop <- psifunC(tau, c(gamprop[1], gam[2]), tau0)
    lanum1 <- sum(lgamma(psi1prop + eta)) - sum(lgamma(psi1prop)) +
      sum(psi1prop)*log(phi) - (gamprop[1] - mugam[1])^2/sig2gam[1]/2
    ladenom1 <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam[1] - mugam[1])^2/sig2gam[1]/2
    probgam1 <- exp(min(lanum1 - ladenom1, 0))
    if(u[1] < probgam1){
      gam[1] <- gamprop[1]
      psi <- psi1prop
      gamaccs[iter,1] <- 1
    }
    psi2prop <- psifunC(tau, c(gam[1], gamprop[2]), tau0)
    lanum2 <- sum(lgamma(psi2prop + eta)) - sum(lgamma(psi2prop)) +
      sum(psi2prop)*log(phi) - (gamprop[2] - mugam[2])^2/sig2gam[2]/2
    ladenom2 <- sum(lgamma(psi + eta)) - sum(lgamma(psi)) +
      sum(psi)*log(phi) - (gam[2] - mugam[2])^2/sig2gam[2]/2
    probgam2 <- exp(min(lanum2 - ladenom2, 0))
    if(u[2] < probgam2){
      gam[2] <- gamprop[2]
      psi <- psi2prop
      gamaccs[iter,2] <- 1
    }
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(nds[k] > 1){
        if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
          lo <- rloggammaC(nds[k], alpha, 1)
          lo <- lo - max(lo)
          ld <- log(D[k]) + lo - log(sum(exp(lo)))
          ldlist[[k]] <- ld
        } else {
          lowlim <- (Dbar[k-1] - D[k-1])/D[k]
          d1 <- rbetarejC(1, alpha, alpha*(nds[k] - 1), lowlim)
          loother <- rloggammaC(nds[k] - 1, alpha, 1)
          ldother <- loother - log(sum(exp(loother)))
          ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
          ldlist[[k]] <- ld
        }
      } else {
        ldlist[[k]] <- log(D[k])
      }
    }
    ldvec <- unlist(ldlist)
    ## ## print("d step finished")
    Nover <- Ntotal - Mtotal
    ## print(c("Nover", Nover))
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, alpha, beta)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, alpha, beta, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), alpha, beta))
      }
      if(Nover > 1){
        ldover <- c(ldover, rloggammaC(Nover - 1, alpha, beta))
      }
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    } else {
      ldover <- rloggammaC(1, alpha, beta)
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    }
    ## ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
       } else {
        lower <- max(M[k+1], N[k]) - N[k]
        psik <- psifunC(tau[k], gam, tau0[k])
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          psik1 <- psifunC(tau[k+1], gam, tau0[k+1])
          etaset <- lower:upper
          etalogprobs <- lgamma(etaset + psik) - lfactorial(etaset) +
            lgamma(upper - etaset + psik1) - lfactorial(upper - etaset)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          if(lower > 0){
            etak <- rktnb(1, psik, lower - 1, psik*(1-phi)/phi)
          } else {
            etak <- rnbinom(1, psik, phi)
          }
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    eta <- diff(N)
    ## ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
        Dbar[k] <- Dbars[min(idxk)]
      } else {
        idxk <- which(Dbars >= D[k])
        Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
      }
    }
    ## alpha step
    alphaprop <- exp(log(alpha) + rnorm(1, 0, alpharwsd))
    u <- runif(1)
    con1 <- Ntotal*log(beta) + Ntotal*mean(ldvec) - balpha
    ## print(c(mean(ldvec), mean(exp(ldvec))))
    laalpha <- (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) -
      Ntotal*(lgamma(alphaprop) - lgamma(alpha))
    probalpha <- exp(min(laalpha, 0))
    u <- runif(1)
    if(u < probalpha){
      alpha <- alphaprop
      alphaaccs[iter,1] <- 1
    }
    ## beta step
    ratebeta <- bbeta + sum(D)
    if(Nover > 0){
      ratebeta <- ratebeta + sum(exp(ldover))
    }
    beta <- rgamma(1, abeta + alpha*Ntotal, ratebeta)
    ## Tuning update step
    if(tune){
      taulogrwsd <- taulogrwsd + rwc*(probtau - accprob[1])
      taurwsd <- exp(taulogrwsd)
      gamlogrwsds <- gamlogrwsds + rwc*(c(probgam1, probgam2) - accprob[2])
      gamrwsds <- exp(gamlogrwsds)
      alphalogrwsd <- alphalogrwsd + rwc*(probalpha - accprob[3])
      alpharwsd <- exp(alphalogrwsd)
    }    
    ## ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alpha, beta, N[-1], tau, Dbar)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alpha = alpha, beta = beta,
                  N = N, tau = tau, Dbar = Dbar)
  rwsds$taulogrwsd <- taulogrwsd
  rwsds$gamlogrwsds <- gamlogrwsds
  rwsds$alphalogrwsd <- alphalogrwsd
  out <- list(draws = draws, initial = initial, rwsds = rwsds, tauaccs = tauaccs,
              gamaccs = gamaccs, alphaaccs = alphaaccs)
  return(out)
}

rbetarej <- function(n, a, b, L){
  if(L > 1){
    warning("lower limit greater than 1")
  }
  if( (a - 1)*(1 - L)^2 + (b - 1)*L^2 >= 0 & b >= 1){
    ## if log concave use ars (if not either a < 1 or b < 1 or both)
    lfun <- function(x, a, b){(a - 1)*log(x) + (b-1)*log(1-x)}
    lfunp <- function(x, a, b){(a-1)/x - (b-1)/(1-x)}
    x <- ars(n=n, f=lfun, fprima=lfunp, x=c(L, L + (1-L)/2, .99),
             xlb=L, xub=1, lb=TRUE, ub=TRUE, a = a, b = b)
  } else {
    prej <- pbeta(L, a, b)
    if(prej >= 0.95 & b < 1) {
      ## if b < 1, use (x-L)(1-L) ~ beta(1, b) proposal if P(rej) is too high
      rej <- 0
      x <- NULL
      lM <- (a < 1)*(a - 1)*log(L)
      while(rej < n){
        y <- rbeta(1, 1, b)
        y <- (1 - L)*y + L
        u <- runif(1)
        if(log(u) < (a - 1)*log(y) - lM){
          rej <- rej + 1
          x <- c(x, y)
        }
      }
    } else if(prej >= 0.95 & a < 1){
      ## if a < 1 and b > 1, use x ~ beta(a, 1)1(x > L) proposal if P(rej) is too high
      rej <- 0
      x <- NULL
      lM <- (b-1)*log(1 - L)
      while(rej < n){
        uy <- runif(1)
        y <- ((1 - L^a)*uy + L^a)^(1/a)
        u <- runif(1)
        if(log(u) < (b-1)*log(1-y) - lM){
          rej <- rej + 1
          x <- c(x, y)
        }
      }
    } else {
      ## if P(rej) isn't too high use naive rejection sampling
      rej <- 0
      x <- NULL
      while(rej < n){
        y <- rbeta(1, a, b)
        if(y > L){
          rej <- rej + 1
          x <- c(x, y)
        }
      }
    }
  }
  return(x)
}

rgammarej <- function(n, alpha, beta, lower){
  rej <- TRUE
  while(rej){
    x <- rgamma(1, alpha, beta)
    if(x > lower){
      rej <- FALSE
    }
  }
  return(x)
}

rloggamma <- function(n, shape, rate){
  log(rgamma(n, shape + 1)) + log(runif(n))/shape - log(rate)
}

psifunR <- function(tau, gam, tau0){
  deltau <- tau - tau0
  if(gam[2] == 0){
    out <- deltau*exp(gam[1])
  } else {
    out <- log(abs(exp(deltau*gam[2]) - 1)) - log(abs(gam[2])) + gam[1] + gam[2]*tau0
    out <- exp(out)
  }
  return(out)
}

logdp <- function(p, con, a){
  out <- con*x - lgamma(x) + (a-1)*log(x)
  return(out)
}

logdpprime <- function(p, con, a){
  out <- con - digamma(x) + (a-1)/x
  return(out)
}

psifunC <- cmpfun(psifunR)
rloggammaC <- cmpfun(rloggamma)
rbetarejC <- cmpfun(rbetarej)
rgammarejC <- cmpfun(rgammarej)
overstaumcmcC <- cmpfun(overstaumcmc)
overftaumcmcC <- cmpfun(overftaumcmc)
staumcmcC <- cmpfun(staumcmc)
ftaumcmcC <- cmpfun(ftaumcmc)
ftautrendmcmcC <- cmpfun(ftautrendmcmc)
overftautrendmcmcC <- cmpfun(overftautrendmcmc)
cgammamcmcC <- cmpfun(cgammamcmc)
