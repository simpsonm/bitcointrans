library(MCMCpack)
library(truncnorm) ## rtruncnorm() is faster than rtnorm() in the msm package
library(ars)
library(compiler)

ftaumcmc <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  nblocks <- length(tau)
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alphad <- initial$alphad
  betad <- initial$betad
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphadlogrwsd <- rwsds$alphadlogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alphadrwsd <- exp(alphadlogrwsd)
  taurejs <- matrix(0, ncol = 1, nrow = niter)
  gamrejs <- matrix(0, ncol = 1, nrow = niter)
  alphadrejs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- rwsds$tune
  lowtarget <- rwsds$lowtarget
  hightarget <- rwsds$hightarget
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alphad, betad
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  lambdanames <- paste("lambda", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alphad", "betad",
                       Nnames, taunames, Dbarnames, lambdanames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    print(iter)
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
    sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    oldsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tau))
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    la <- (oldsqdiff - newsqdiff)/(2*sigma2) +
        sum((deltaprop - delta)*(gam*log(phi*lambda) - 1/10)) +
        sum(lgamma(gam*delta) - lgamma(gam*deltaprop) + zetaprop - zetaold)
    u <- runif(1)
    if(log(u) < la){
      delta <- deltaprop
      tau <- tauprop
    } else {
      taurejs[iter] <- 1
    }
    ## lambda step
    eta <- diff(N)
    lambda <- rgamma(nblocks, eta + gam*delta, phi + 1)
    ## phi step
    phi <- rgamma(1, aphi + gam*sum(delta), bphi + sum(lambda))
    ## print("lambda and phi steps finished")
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    if(gamprop > 0){
      u <- runif(1)
      con1 <- (-bgam + sum(delta*log(phi*lambda)))
      lanum <-   gamprop*con1 + (agam - 1)*log(gamprop) - sum(lgamma(gamprop*delta))
      ladenom <- gam*con1     + (agam - 1)*log(gam)     - sum(lgamma(gam*delta))
      lagam <- lanum - ladenom
      if(log(u) < lagam){
        gam <- gamprop
      } else {
        gamrejs[iter,1] <- 1
      }
    } else {
      gamrejs[iter,1] <- 1
    }
    ## print("gamma step finished")
    ## print("Dbar step finished")
    ## alphad step - Metrop
    alphadprop <- rnorm(1, alphad, alphadrwsd)
    if(alphadprop > 0){
      u <- runif(1)
      con1 <- -balpha + Ntotal*log(betad) + sum(eta*log(D))
      laalpha <- (alphadprop - alphad)*con1 + (aalpha - 1)*log(alphadprop/alphad) -
        sum(lgamma(eta*alphadprop) - lgamma(eta*alphad))
      if(log(u) < laalpha){
        alphad <- alphadprop
      } else {
        alphadrejs[iter,1] <- 1
      }
    } else {
      alphadrejs[iter,1] <- 1
    }
    ## print("alphad step finished")
    ## betad step
    betad <- rgamma(1, abeta + alphad*Ntotal, bbeta + sum(D))
    ## Tuning update step
    ## print("betad step finished")
    if(tune & iter %% H == 0){
      tauacc <- 1 - mean(taurejs[(iter - H + 1):iter,])
      taulogrwsd <- taulogrwsd + ((tauacc > hightarget[1]) - (tauacc < lowtarget[1]))*rwc
      taurwsd <- exp(taulogrwsd)
      gamacc <- 1 - mean(gamrejs[(iter - H + 1):iter,])
      gamlogrwsd <- gamlogrwsd + ((gamacc > hightarget[2]) - (gamacc < lowtarget[2]))*rwc
      gamrwsd <- exp(gamlogrwsd)
      alphadacc <- 1 - mean(alphadrejs[(iter - H + 1):iter,])
      alphadlogrwsd <- alphadlogrwsd + ((alphadacc > hightarget[3]) - (alphadacc < lowtarget[3]))*rwc
      alphadrwsd <- exp(alphadlogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alphad, betad, N[-1], tau, Dbar, lambda)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alphad = alphad, betad = betad,
                  N = N, tau = tau, Dbar = Dbar, lambda = lambda)
  rwsds <- list(taulogrwsd = taulogrwsd, gamlogrwsd = gamlogrwsd, alphadlogrwsd = alphadlogrwsd)
  out <- list(draws = draws, initial = initial, rwsds = rwsds, taurejs = taurejs,
              gamrejs = gamrejs, alphadrejs = alphadrejs)
  return(out)
}

overftaumcmc <- function(niter, data, initial, prior, rwsds, tune, rwc, H, tauchol, hightarget, lowtarget){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  nblocks <- length(tau)
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alphad <- initial$alphad
  betad <- initial$betad
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsd <- rwsds$gamlogrwsd
  alphadlogrwsd <- rwsds$alphadlogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsd <- exp(gamlogrwsd)
  alphadrwsd <- exp(alphadlogrwsd)
  taurejs <- matrix(0, ncol = 1, nrow = niter)
  gamrejs <- matrix(0, ncol = 1, nrow = niter)
  alphadrejs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- rwsds$tune
  lowtarget <- rwsds$lowtarget
  hightarget <- rwsds$hightarget
  tauchol <- rwsds$tauchol
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alphad, betad
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  lambdanames <- paste("lambda", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alphad", "betad",
                       Nnames, taunames, Dbarnames, lambdanames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    print(iter)
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
    sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## print("sigma step finished")
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    oldsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tau))
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    la <- (oldsqdiff - newsqdiff)/(2*sigma2) +
        sum((deltaprop - delta)*(gam*log(phi*lambda) - 1/10)) +
        sum(lgamma(gam*delta) - lgamma(gam*deltaprop) + zetaprop - zetaold)
    u <- runif(1)
    if(log(u) < la){
      delta <- deltaprop
      tau <- tauprop
    } else {
      taurejs[iter] <- 1
    }
    ## lambda step
    eta <- diff(N)
    lambda <- rgamma(nblocks, eta + gam*delta, phi + 1)
    ## phi step
    phi <- rgamma(1, aphi + gam*sum(delta), bphi + sum(lambda))
    ## print("lambda and phi steps finished")
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    if(gamprop > 0){
      u <- runif(1)
      con1 <- (-bgam + sum(delta*log(phi*lambda)))
      lanum <-   gamprop*con1 + (agam - 1)*log(gamprop) - sum(lgamma(gamprop*delta))
      ladenom <- gam*con1     + (agam - 1)*log(gam)     - sum(lgamma(gam*delta))
      lagam <- lanum - ladenom
      if(log(u) < lagam){
        gam <- gamprop
      } else {
        gamrejs[iter,1] <- 1
      }
    } else {
      gamrejs[iter,1] <- 1
    }
    ## print("gamma step finished")
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
        lo <- rloggammaC(nds[k], alphad, 1)
        ld <- log(D[k]) + lo - log(sum(exp(lo)))
        ldlist[[k]] <- ld
      } else {
        lowlim <- (Dbar[k-1] - D[k-1])/D[k]
        ## print("rej step")
        ## print(k)
        d1 <- rbetarejC(1, alphad, alphad*(nds[k] - 1), lowlim)
        ## print("rej step finished")
        loother <- rloggammaC(nds[k] - 1, alphad, 1)
        ldother <- loother - log(sum(exp(loother)))
        ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
        ldlist[[k]] <- ld
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    Nover <- Ntotal - Mtotal
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, alphad, betad)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, alphad, betad, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), alphad, betad))
      }
      if(Nover > 1){
        ldover <- c(ldover, rloggammaC(Nover - 1, alphad, betad))
      }
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    } else {
      ldlist[[nblocks + 1]] <- -Inf
    }
    ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
      } else {
        lower <- M[k+1] - N[k]
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          etaset <- lower:upper
          etalogprobs <- dpois(etaset, lambda[k], log = TRUE)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          etak <- rpoistrunlowC(1, lambda[k], lower)
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
      } else {
        idxk <- which(Dbars >= D[k])
      }
      Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
    }
    ## print("Dbar step finished")
    ## alphad betad step - Metrop
    alphadprop <- rnorm(1, alphad, alphadrwsd)
    if(alphadprop > 0){
      u <- runif(1)
      con1 <- Ntotal*mean(ldvec) - Ntotal*log(bbeta + Ntotal*mean(exp(ldvec))) - balpha
      laalpha <- (alphadprop - alphad)*con1 + (aalpha - 1)*log(alphadprop/alphad) -
        Ntotal*(lgamma(alphadprop) - lgamma(alphad)) +
        lgamma(alphadprop*Ntotal + abeta - 1) - lgamma(alphad*Ntotal + abeta - 1)
      if(log(u) < laalpha){
        alphad <- alphadprop
        betad <- rgamma(1, abeta + alphad*Ntotal, bbeta + Ntotal*mean(exp(ldvec)))
      } else {
        alphadrejs[iter,1] <- 1
      }
    } else {
      alphadrejs[iter,1] <- 1
    }
    ## Tuning update step
    if(tune & iter %% H == 0){
      tauacc <- 1 - mean(taurejs[(iter - H + 1):iter,])
      taulogrwsd <- taulogrwsd + ((tauacc > hightarget[1]) - (tauacc < lowtarget[1]))*rwc
      taurwsd <- exp(taulogrwsd)
      gamacc <- 1 - mean(gamrejs[(iter - H + 1):iter,])
      gamlogrwsd <- gamlogrwsd + ((gamacc > hightarget[2]) - (gamacc < lowtarget[2]))*rwc
      gamrwsd <- exp(gamlogrwsd)
      alphadacc <- 1 - mean(alphadrejs[(iter - H + 1):iter,])
      alphadlogrwsd <- alphadlogrwsd + ((alphadacc > hightarget[3]) - (alphadacc < lowtarget[3]))*rwc
      alphadrwsd <- exp(alphadlogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alphad, betad, N[-1], tau, Dbar, lambda)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alphad = alphad, betad = betad,
                  N = N, tau = tau, Dbar = Dbar, lambda = lambda)
  rwsds <- list(taulogrwsd = taulogrwsd, gamlogrwsd = gamlogrwsd, alphadlogrwsd = alphadlogrwsd)
  out <- list(draws = draws, initial = initial, rwsds = rwsds, taurejs = taurejs,
              gamrejs = gamrejs, alphadrejs = alphadrejs)
  return(out)
}

overstaumcmc <- function(niter, data, initial, prior, rwsds, tune, rwc, H, hightarget, lowtarget){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  nblocks <- length(tau)
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alphad <- initial$alphad
  betad <- initial$betad
  ## RW sds and bookkeeping
  taulogrwsds <- rwsds$taulogrwsds
  gamlogrwsd <- rwsds$gamlogrwsd
  alphadlogrwsd <- rwsds$alphadlogrwsd
  taurwsds <- exp(taulogrwsds)
  gamrwsd <- exp(gamlogrwsd)
  alphadrwsd <- exp(alphadlogrwsd)
  taurejs <- matrix(0, ncol = nblocks, nrow = niter)
  gamrejs <- matrix(0, ncol = 1, nrow = niter)
  alphadrejs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- rwsds$tune
  lowtarget <- rwsds$lowtarget
  hightarget <- rwsds$hightarget
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alphad, betad
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  lambdanames <- paste("lambda", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alphad", "betad",
                       Nnames, taunames, Dbarnames, lambdanames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    print(iter)
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
    sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
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
      oldsqdiff <- mapply(function(x, y){sum((x - y)^2)}, tfull, tau)
      newsqdiff <- mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop)
      la <- (sum(oldsqdiff) - sum(newsqdiff))/(2*sigma2) +
        (deltakprop - deltakold)*(gam*log(phi*lambda[k]) - 1/10) +
        lgamma(gam*deltakold) - lgamma(gam*deltakprop) + zetakprop - zetakold
      u <- runif(1)
      if(log(u) < la){
        delta[k] <- deltakprop
        tau <- tauprop
      } else {
        taurejs[iter, k] <- 1
      }
    }
    tau0 <- c(0, tau)
    delta <- diff(tau0)
    ## lambda step
    eta <- diff(N)
    lambda <- rgamma(nblocks, eta + gam*delta, phi + 1)
    ## phi step
    phi <- rgamma(1, aphi + gam*sum(delta), bphi + sum(lambda))
    ## print("lambda and phi steps finished")
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    if(gamprop > 0){
      u <- runif(1)
      con1 <- (-bgam + sum(delta*log(phi*lambda)))
      lanum <-   gamprop*con1 + (agam - 1)*log(gamprop) - sum(lgamma(gamprop*delta))
      ladenom <- gam*con1     + (agam - 1)*log(gam)     - sum(lgamma(gam*delta))
      lagam <- lanum - ladenom
      if(log(u) < lagam){
        gam <- gamprop
      } else {
        gamrejs[iter,1] <- 1
      }
    } else {
      gamrejs[iter,1] <- 1
    }
    ## print("gamma step finished")
    ## d step
    ldlist <- list()
    for(k in 1:nblocks){
      if(N[k] == M[k]){  ## note: N[k] == N_{k+1} and M[k] == M_{k+1}
        lo <- rloggammaC(nds[k], alphad, 1)
        ldlist[[k]] <- lo + log(D[k]) - log(sum(exp(lo)))
      } else {
        lowlim <- (Dbar[k-1] - D[k-1])/D[k]
        d1 <- rbetarejC(1, alphad, alphad*(nds[k] - 1), lowlim)
        loother <- rloggammaC(nds[k] - 1, alphad, 1)
        ldother <- loother - log(sum(exp(loother)))
        ld <- c(log(d1), log((1 - d1)) + ldother) + log(D[k])
        ldlist[[k]] <- ld
      }
    }
    ldvec <- unlist(ldlist)
    ## print("d step finished")
    Nover <- Ntotal - Mtotal
    if(Nover > 0){
      lowlim <- Dbar[nblocks] - D[nblocks]
      prej <- pgamma(lowlim, alphad, betad)
      if(prej < 0.95){
        ldover <- log(rgammarejC(1, alphad, betad, lowlim))
      } else{
        u <- runif(1)
        ldover <- log(qgamma(prej + u * (1 - prej), alphad, betad))
      }
      if(Nover > 1){
        ldover <- c(ldover, rloggammaC(Nover-1, alphad, betad))
      }
      ldlist[[nblocks + 1]] <- ldover
      ldvec <- c(ldvec, ldover)
    } else {
      ldlist[[nblocks + 1]] <- Inf
    }
    ## print("Nover d step finished")
    ## N step
    for(k in 1:nblocks){
      if(D[k] + exp(ldlist[[k+1]][1]) <= Dbar[k]){
        N[k+1] = M[k+1]
      } else {
        lower <- M[k+1] - N[k]
        if(k < nblocks){
          upper <- N[k + 2] - N[k]
          etaset <- lower:upper
          etalogprobs <- dpois(etaset, lambda[k], log = TRUE)
          etaprobs <- exp(etalogprobs - max(etalogprobs))
          etak <- sample(etaset, 1, prob = etaprobs)
        } else {
          etak <- rpoistrunlowC(1, lambda[k], lower)
        }
        N[k+1] <- N[k] + etak
      }
    }
    Ntotal <- N[nblocks + 1]
    ## print("N step finished")
    ## Dbar step
    for(k in 1:nblocks){
      if(N[k+1] > M[k+1]){
        idxk <- which(Dbars >= D[k] & Dbars <= D[k] + exp(ldlist[[k+1]][1]))
      } else {
        idxk <- which(Dbars >= D[k])
      }
      Dbar[k] <- Dbars[idxk[sample.int(length(idxk), 1)]]
    }
    ## alphad betad step - Metrop
    alphadprop <- rnorm(1, alphad, alphadrwsd)
    if(alphadprop > 0){
      u <- runif(1)
      con1 <- Ntotal*mean(ldvec) - Ntotal*log(bbeta + Ntotal*mean(exp(ldvec))) - balpha
      laalpha <- (alphadprop - alphad)*con1 + (aalpha - 1)*log(alphadprop/alphad) -
        Ntotal*(lgamma(alphadprop) - lgamma(alphad)) +
        lgamma(alphadprop*Ntotal + abeta - 1) - lgamma(alphad*Ntotal + abeta - 1)
      if(log(u) < laalpha){
        alphad <- alphadprop
        betad <- rgamma(1, abeta + alphad*Ntotal, bbeta + Ntotal*mean(exp(ldvec)))
      } else {
        alphadrejs[iter,1] <- 1
      }
    } else {
      alphadrejs[iter,1] <- 1
    }
    ## Tuning update step
    ## print("betad step finished")
    if(tune & iter %% H == 0){
      tauacc <- 1 - apply(taurejs[(iter - H + 1):iter,], 2, mean)
      taulogrwsds <- taulogrwsds + ((tauacc > hightarget[1]) - (tauacc < lowtarget[1]))*rwc
      taurwsds <- exp(taulogrwsds)
      gamacc <- 1 - mean(gamrejs[(iter - H + 1):iter,])
      gamlogrwsd <- gamlogrwsd + ((gamacc > hightarget[2]) - (gamacc < lowtarget[2]))*rwc
      gamrwsd <- exp(gamlogrwsd)
      alphadacc <- 1 - mean(alphadrejs[(iter - H + 1):iter,])
      alphadlogrwsd <- alphadlogrwsd + ((alphadacc > hightarget[3]) - (alphadacc < lowtarget[3]))*rwc
      alphadrwsd <- exp(alphadlogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alphad, betad, N[-1], tau, Dbar, lambda)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alphad = alphad, betad = betad,
                  N = N, tau = tau, Dbar = Dbar, lambda = lambda)
  rwsds <- list(taulogrwsds = taulogrwsds, gamlogrwsd = gamlogrwsd, alphadlogrwsd = alphadlogrwsd)
  out <- list(draws = draws, initial = initial, rwsds = rwsds, taurejs = taurejs,
              gamrejs = gamrejs, alphadrejs = alphadrejs)
  return(out)
}

staumcmc <- function(niter, data, initial, prior, rwsds, tune, rwc, H, hightarget, lowtarget){
  ## Initial values
  tau <- initial$tau
  tau0 <- c(0, tau)
  delta <- diff(tau0)
  nblocks <- length(tau)
  sigma <- initial$sigma
  lambda <- initial$lambda
  gam <- initial$gam
  phi <- initial$phi
  N <- initial$N
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alphad <- initial$alphad
  betad <- initial$betad
  ## RW sds and bookkeeping
  taulogrwsds <- rwsds$taulogrwsds
  gamlogrwsd <- rwsds$gamlogrwsd
  alphadlogrwsd <- rwsds$alphadlogrwsd
  taurwsds <- exp(taulogrwsds)
  gamrwsd <- exp(gamlogrwsd)
  alphadrwsd <- exp(alphadlogrwsd)
  taurejs <- matrix(0, ncol = nblocks, nrow = niter)
  gamrejs <- matrix(0, ncol = 1, nrow = niter)
  alphadrejs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- rwsds$tune
  lowtarget <- rwsds$lowtarget
  hightarget <- rwsds$hightarget
  ## Priors
  vsigma <- prior$vsigma
  ssigma <- prior$ssigma
  aphi <- prior$aphi
  bphi <- prior$bphi
  agam <- prior$agamma
  bgam <- prior$bgamma
  Dbars <- prior$Dbars
  aalpha <- prior$aalpha
  balpha <- prior$balpha
  abeta <- prior$abeta
  bbeta <- prior$bbeta
  aomega <- (vsigma + 1)/2
  bomega <- 1/(vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alphad, betad
  draws <- matrix(0, nrow = niter, ncol = nblocks + 1 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  lambdanames <- paste("lambda", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", "gam", "alphad", "betad",
                       Nnames, taunames, Dbarnames, lambdanames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    print(iter)
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
    sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
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
      oldsqdiff <- mapply(function(x, y){sum((x - y)^2)}, tfull, tau)
      newsqdiff <- mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop)
      la <- (sum(oldsqdiff) - sum(newsqdiff))/(2*sigma2) +
        (deltakprop - deltakold)*(gam*log(phi*lambda[k]) - 1/10) +
        lgamma(gam*deltakold) - lgamma(gam*deltakprop) + zetakprop - zetakold
      u <- runif(1)
      if(log(u) < la){
        delta[k] <- deltakprop
        tau <- tauprop
      } else {
        taurejs[iter, k] <- 1
      }
    }
    tau0 <- c(0, tau)
    delta <- diff(tau0)
    ## lambda step
    eta <- diff(N)
    lambda <- rgamma(nblocks, eta + gam*delta, phi + 1)
    ## phi step
    phi <- rgamma(1, aphi + gam*sum(delta), bphi + sum(lambda))
    ## print("lambda and phi steps finished")
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(1, gam, gamrwsd)
    if(gamprop > 0){
      u <- runif(1)
      con1 <- (-bgam + sum(delta*log(phi*lambda)))
      lanum <-   gamprop*con1 + (agam - 1)*log(gamprop) - sum(lgamma(gamprop*delta))
      ladenom <- gam*con1     + (agam - 1)*log(gam)     - sum(lgamma(gam*delta))
      lagam <- lanum - ladenom
      if(log(u) < lagam){
        gam <- gamprop
      } else {
        gamrejs[iter,1] <- 1
      }
    } else {
      gamrejs[iter,1] <- 1
    }
    ## print("gamma step finished")
    ## alphad step - Metrop
    alphadprop <- rnorm(1, alphad, alphadrwsd)
    if(alphadprop > 0){
      u <- runif(1)
      con1 <- -balpha + Ntotal*log(betad) + sum(eta*log(D))
      laalpha <- (alphadprop - alphad)*con1 + (aalpha - 1)*log(alphadprop/alphad) -
        sum(lgamma(eta*alphadprop) - lgamma(eta*alphad))
      if(log(u) < laalpha){
        alphad <- alphadprop
      } else {
        alphadrejs[iter,1] <- 1
      }
    } else {
      alphadrejs[iter,1] <- 1
    }
    ## print("alphad step finished")
    ## betad step
    betad <- rgamma(1, abeta + alphad*Ntotal, bbeta + sum(D))
    ## Tuning update step
    ## print("betad step finished")
    if(tune & iter %% H == 0){
      tauacc <- 1 - apply(taurejs[(iter - H + 1):iter,], 2, mean)
      taulogrwsds <- taulogrwsds + ((tauacc > hightarget[1]) - (tauacc < lowtarget[1]))*rwc
      taurwsds <- exp(taulogrwsds)
      gamacc <- 1 - mean(gamrejs[(iter - H + 1):iter,])
      gamlogrwsd <- gamlogrwsd + ((gamacc > hightarget[2]) - (gamacc < lowtarget[2]))*rwc
      gamrwsd <- exp(gamlogrwsd)
      alphadacc <- 1 - mean(alphadrejs[(iter - H + 1):iter,])
      alphadlogrwsd <- alphadlogrwsd + ((alphadacc > hightarget[3]) - (alphadacc < lowtarget[3]))*rwc
      alphadrwsd <- exp(alphadlogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alphad, betad, N[-1], tau, Dbar, lambda)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alphad = alphad, betad = betad,
                  N = N, tau = tau, Dbar = Dbar, lambda = lambda)
  rwsds <- list(taulogrwsds = taulogrwsds, gamlogrwsd = gamlogrwsd, alphadlogrwsd = alphadlogrwsd)
  out <- list(draws = draws, initial = initial, rwsds = rwsds, taurejs = taurejs,
              gamrejs = gamrejs, alphadrejs = alphadrejs)
  return(out)
}

rbetarej <- function(n, a, b, L){
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

rpoistrunlow <- function(n, lambda, lower, nup = 200){
  lprobs <- dpois(lower:max(ceiling(lambda), lower + nup), lambda, TRUE)
  mlprobs <- max(lprobs)
  k <- 0
  while(lprobs[nup + k] - mlprobs > -700){
    k <- k+1
    lprobs <- c(lprobs, dpois(lower + k, lambda, TRUE))
  }
  probs <- exp(lprobs - max(lprobs))
  out <- sample(length(probs), n, prob = probs)
  return(out + lower - 1)
}

rloggamma <- function(n, shape, rate){
  log(rgamma(n, shape + 1)) + log(runif(n))/shape - log(rate)
}

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
  Ntotal <- N[nblocks + 1]
  Dbar <- initial$Dbar
  alphad <- initial$alphad
  betad <- initial$betad
  psi <- psifun(tau, gam, tau0)
  ## RW sds and bookkeeping
  taulogrwsd <- rwsds$taulogrwsd
  gamlogrwsds <- rwsds$gamlogrwsds
  alphadlogrwsd <- rwsds$alphadlogrwsd
  taurwsd <- exp(taulogrwsd)
  gamrwsds <- exp(gamlogrwsds)
  alphadrwsd <- exp(alphadlogrwsd)
  taurejs <- matrix(0, ncol = 1, nrow = niter)
  gamrejs <- matrix(0, ncol = 2, nrow = niter)
  alphadrejs <- matrix(0, ncol = 1, nrow = niter)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- rwsds$tune
  lowtarget <- rwsds$lowtarget
  hightarget <- rwsds$hightarget
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
  bomega <- 1/(vsigma*ssigma^2)
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  Mtotal <- M[nblocks + 1]
  D <- data$D
  nds <- diff(M) ## note: M[1] = 0 and N[1] = 0, M_k = M[k + 1]
  ## Store draws: tau, sigma, lambda, gam, phi, N, Dbar, alphad, betad
  draws <- matrix(0, nrow = niter, ncol = nblocks + 2 + nblocks + 1 + 1 + nblocks + nblocks + 2)
  taunames <- paste("tau", 1:nblocks, sep = "")
  lambdanames <- paste("lambda", 1:nblocks, sep = "")
  Nnames <- paste("N", 1:nblocks, sep = "")
  Dbarnames <- paste("Dbar", 1:nblocks, sep = "")
  colnames(draws) <- c("sigma", "phi", paste("gam", 1:2, sep=""), "alphad", "betad",
                       Nnames, taunames, Dbarnames, lambdanames)
  ## Initialize DA stuff
  nrej <- rep(0, nblocks)
  tfull <- list()
  ## Main loop
  for(iter in 1:niter){
    print(iter)
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
    sigma2 <- 1/rgamma(1, asigma, bsigma/2 + omega )
    sigma <- sqrt(sigma2)
    ## tau step
    zetaold <- log(delta)
    zetaprop <- zetaold + taurwsd*tauchol%*%rnorm(nblocks)
    deltaprop <- exp(zetaprop)
    tauprop <- cumsum(deltaprop)
    tau0prop <- c(0, tauprop[-nblocks])
    oldsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tau))
    newsqdiff <- sum(mapply(function(x, y){sum((x - y)^2)}, tfull, tauprop))
    psiprop <- psifunC(tauprop, gam, tau0prop)
    la <- (oldsqdiff - newsqdiff)/(2*sigma2) +
      sum((psiprop - psi)*log(phi*lambda)) + (tau[nblocks] - tauprop[nblocks])/10 +
      sum(lgamma(psi)) - sum(lgamma(psiprop)) + sum(zetaprop) - sum(zetaold)
    u <- runif(1)
    if(log(u) < la){
      delta <- deltaprop
      tau <- tauprop
      psi <- psiprop
      tau0 <- tau0prop
    } else {
      taurejs[iter] <- 1
    }
    ## lambda step
    eta <- diff(N)
    lambda <- rgamma(nblocks, eta + psi, phi + 1)
    ## phi step
    phi <- rgamma(1, aphi + sum(psi), bphi + sum(lambda))
    ## print("lambda and phi steps finished")
    ## gamma step
    ## RW Metropolis with tuned proposals
    gamprop <- rnorm(2, gam, gamrwsds)
    u <- runif(2)
    psi1prop <- psifunC(tau, c(gamprop[1], gam[2]), tau0)
    print("psi1prop")
    print(summary(psi1prop))
    lanum1 <- sum(psi1prop*log(phi*lambda)) - sum(lgamma(psi1prop)) +
      (gamprop[1] - mugam[1])^2/sig2gam[1]
    ladenom1 <- sum(psi*log(phi*lambda)) - sum(lgamma(psi)) +
      (gam[1] - mugam[1])^2/sig2gam[1]
    lagam1 <- lanum1 - ladenom1
    if(log(u[1]) < lagam1){
      gam[1] <- gamprop[1]
      psi <- psi1prop
    } else {
      gamrejs[iter,1] <- 1
    }
    psi2prop <- psifunC(tau, c(gam[1], gamprop[2]), tau0)
    print("psi2prop")
    print(summary(psi2prop))
    lanum2 <- sum(psi2prop*log(phi*lambda)) - sum(lgamma(psi2prop)) +
      (gamprop[2] - mugam[2])^2/sig2gam[2]
    ladenom2 <- sum(psi*log(phi*lambda)) - sum(lgamma(psi)) +
      (gam[2] - mugam[2])^2/sig2gam[2]
    lagam2 <- lanum2 - ladenom2
    if(log(u[2]) < lagam2){
      gam[2] <- gamprop[2]
      psi <- psi2prop
    } else {
      gamrejs[iter,2] <- 1
    }
    ## alphad step - Metrop
    alphadprop <- rnorm(1, alphad, alphadrwsd)
    if(alphadprop > 0){
      u <- runif(1)
      con1 <- -balpha + Ntotal*log(betad) + sum(eta*log(D))
      laalpha <- (alphadprop - alphad)*con1 + (aalpha - 1)*log(alphadprop/alphad) -
        sum(lgamma(eta*alphadprop) - lgamma(eta*alphad))
      if(log(u) < laalpha){
        alphad <- alphadprop
      } else {
        alphadrejs[iter,1] <- 1
      }
    } else {
      alphadrejs[iter,1] <- 1
    }
    ## print("alphad step finished")
    ## betad step
    betad <- rgamma(1, abeta + alphad*Ntotal, bbeta + sum(D))
    ## Tuning update step
    ## print("betad step finished")
    if(tune & iter %% H == 0){
      tauacc <- 1 - mean(taurejs[(iter - H + 1):iter,])
      taulogrwsd <- taulogrwsd + ((tauacc > hightarget[1]) - (tauacc < lowtarget[1]))*rwc
      taurwsd <- exp(taulogrwsd)
      gamacc <- 1 - apply(gamrejs[(iter - H + 1):iter,], 2, mean)
      gamlogrwsds <- gamlogrwsds + ((gamacc > hightarget[2]) - (gamacc < lowtarget[2]))*rwc
      gamrwsds <- exp(gamlogrwsds)
      alphadacc <- 1 - mean(alphadrejs[(iter - H + 1):iter,])
      alphadlogrwsd <- alphadlogrwsd + ((alphadacc > hightarget[3]) - (alphadacc < lowtarget[3]))*rwc
      alphadrwsd <- exp(alphadlogrwsd)
    }
    ## print("tune step finished")
    draws[iter,] <- c(sigma, phi, gam, alphad, betad, N[-1], tau, Dbar, lambda)
  }
  ## Collect draws and everything else for output
  initial <- list(sigma = sigma, phi = phi, gam = gam, alphad = alphad, betad = betad,
                  N = N, tau = tau, Dbar = Dbar, lambda = lambda)
  rwsds <- list(taulogrwsd = taulogrwsd, gamlogrwsds = gamlogrwsds, alphadlogrwsd = alphadlogrwsd)
  out <- list(draws = draws, initial = initial, rwsds = rwsds, taurejs = taurejs,
              gamrejs = gamrejs, alphadrejs = alphadrejs)
  return(out)
}

psifun <- function(tau, gam, tau0){
  if(gam[2] == 0){
    out <- (tau - tau0)*exp(gam[1])
  } else {
    out <- (exp(gam[1] + tau*gam[2]) - exp(gam[1] + tau0*gam[2]))/gam[2]
  }
  return(out)
}

psifunC <- cmpfun(psifun)
rloggammaC <- cmpfun(rloggamma)
rpoistrunlowC <- cmpfun(rpoistrunlow)
rbetarejC <- cmpfun(rbetarej)
rgammarejC <- cmpfun(rgammarej)
overstaumcmcC <- cmpfun(overstaumcmc)
overftaumcmcC <- cmpfun(overftaumcmc)
staumcmcC <- cmpfun(staumcmc)
ftaumcmcC <- cmpfun(ftaumcmc)
ftautrendmcmcC <- cmpfun(ftautrendmcmc)
