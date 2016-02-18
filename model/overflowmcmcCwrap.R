Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
## required so that armadillo uses c++11 for gamma random number generation
library(MCMCpack)
library(truncnorm) ## rtruncnorm() is faster than rtnorm() in the msm package
library(ars)
library(compiler)
library(aster)
library(Rcpp)
library(RcppArmadillo)
source("overflowmcmc.R")
sourceCpp("overflowmcmcC.cpp")

overtrendmcmcCwrap <- function(niter, data, initial, prior, rwsds){
  ## Initial values
  tau <- as.numeric(c(initial$tau))
  nblocks <- length(tau)
  sigma <- as.numeric(initial$sigma)
  gam <- as.numeric(c(initial$gam))
  phi <- as.numeric(initial$phi)
  N <- as.numeric(c(initial$N))
  Dbar <- as.numeric(c(initial$Dbar))
  alpha <- as.numeric(initial$alpha)
  beta <- as.numeric(initial$beta)
  ## RW sds and bookkeeping
  taulogrwsd <- as.numeric(rwsds$taulogrwsd)
  gamlogrwsds <- as.numeric(rwsds$gamlogrwsds)
  alphalogrwsd <- as.numeric(rwsds$alphalogrwsd)
  rwc <- rwsds$rwc
  H <- rwsds$H
  tune <- (rwsds$tune)*1
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
  ## Data
  tobs <- data$tobs
  m <- data$m
  M <- data$M
  D <- data$D
  out <- overtrendmcmcCPP(niter, tobs, D, M, m, tau, sigma, gam, phi, N, Dbar, alpha, beta,
                        vsigma, ssigma, aphi, bphi, mugam, sig2gam, Dbars, aalpha, balpha,
                        abeta, bbeta, taulogrwsd, gamlogrwsds, alphalogrwsd,
                        rwc, H, tune, lowtarget, hightarget, tauchol, rtruncnorm, rktnb,
                        rbetarejC)
  return(out)
}

## List overtrendmcmcC(int niter, arma::vec tobs, arma::vec D, arma::vec M, arma::vec m,
## 			 arma::vec tau, double sigma, arma::vec gam, double phi, arma::vec N,
## 			 arma::vec Dbar, double alpha, double beta, 
## 			 double vsigma, double ssigma, double aphi, double bphi, 
## 			 arma::vec mugam, arma::vec sig2gam, arma::vec Dbars, double aalpha,
## 			 double balpha, double abeta, double bbeta,
## 			 double taulogrwsd, arma::vec gamlogrwsds, double alphalogrwsd, 
## 			 double rwc, int H, int tune,
## 			 arma::vec lowtarget, arma::vec hightarget, arma::mat tauchol,
## 			 Function rtruncnorm, Function rktnb, Function rbetarejC){

