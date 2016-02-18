/* To do:
   4) several functions called from R - need to written in C++
   5) N_K = M_K must be enforced. Needs to be fixed.
*/

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;

// used for generating log gamma distributed random variates when shape < 1
arma::vec rloggamma(int n, double shape, double scale){
    arma::vec temp = log(arma::randg(n, arma::distr_param(shape + 1.0, 1.0))) + 
      log(arma::randu(n))/shape + log(scale);
  return temp;
}

// simple rejection sampler for gamma distributed random variates
arma::vec rgammarej(double shape, double scale, double lower){
  arma::vec x = arma::randg(1, arma::distr_param(shape, scale));
  while(x(0) < lower){
    x = arma::randg(1, arma::distr_param(shape, scale));
  }
  return x;
}

// compute psi vector
// [[Rcpp::export]]
arma::vec psifun(arma::vec tau, arma::vec gam, arma::vec tau0){
  arma::vec deltau = tau - tau0;
  arma::vec out(tau.size());
  if(gam(1) == 0){
    out = deltau*exp(gam(0));
  } else {
    // abs(arma::vec) works as expected
    // abs(double) coerces to integer, causing lots of problems
    // std::abs(double)
    out = log(abs(exp(deltau*gam(1)) - 1)) + gam(0) + gam(1)*tau0 - log(std::abs(gam(1)));
    out = exp(out);
  }
  return out;
}

// [[Rcpp::export]]
List overtrendmcmcCPP(int niter, arma::vec tobs, arma::vec D, arma::uvec M, arma::vec m,
			 arma::vec tau, double sigma, arma::vec gam, double phi, arma::uvec N,
			 arma::vec Dbar, double alpha, double beta, 
			 double vsigma, double ssigma, double aphi, double bphi, 
			 arma::vec mugam, arma::vec sig2gam, arma::vec Dbars, double aalpha,
			 double balpha, double abeta, double bbeta,
			 double taulogrwsd, arma::vec gamlogrwsds, double alphalogrwsd, 
			 double rwc, int H, int tune,
			 arma::vec lowtarget, arma::vec hightarget, arma::mat tauchol,
			 Function rtruncnorm, Function rktnb, Function rbetarejC){
  Rprintf("Initializing 1");
  int nblocks = tau.size();
  arma::vec tau0full(nblocks + 1); tau0full(0) = 0; tau0full.subvec(1, nblocks) = tau;
  arma::vec tau0 = tau0full.subvec(0, nblocks - 1);
  arma::vec delta = diff(tau0full);
  arma::uvec eta = diff(N);
  int etak;
  arma::vec psi = psifun(tau, gam, tau0);
  double taurwsd = exp(taulogrwsd);
  arma::vec gamrwsds = exp(gamlogrwsds);
  double alpharwsd = exp(alphalogrwsd);
  arma::mat taurejs(niter, 1);
  arma::mat gamrejs(niter, 2);
  arma::mat alpharejs(niter, 1);
  Rprintf(" 2");  
  double aomega = (vsigma + 1)/2;
  double bomega = 1/(vsigma*pow(ssigma, 2.0));
  double asigma;
  double con2 = aalpha - 1;
  arma::uvec nds = diff(M);  //////////////////////
  arma::mat draws(niter, 3*nblocks + 6);
  arma::uvec nrej(nblocks);
  arma::vec rho(nblocks);
  arma::uvec nrejcumsum(nblocks);
  Rprintf(" 3"); 
  double bsigma;
  double omega;
  double sigma2;
  arma::vec zeta = log(delta);
  arma::vec zetaprop(nblocks);
  arma::vec deltaprop(nblocks);
  arma::vec tauprop(nblocks);
  arma::vec tau0prop(nblocks); tau0prop(0) = 0;
  arma::vec psiprop(nblocks);
  double lataunum;
  double lataudenom;
  double gamprop;
  arma::vec gamtemp(2);
  double lanum;
  Rprintf(" 4"); 
  double ladenom;
  int Nover;
  double lowlim;
  double prej;
  double d1;
  int lower;
  int upper;
  arma::vec dk1(nblocks);
  Rprintf(" 5\n");
  double mnldvec;
  double mnexpldvec;
  double alphaprop;
  int idxklow;
  int idxkhigh;
  double scalebeta = 1/(bbeta + sum(D));
  Rprintf("Beginning main loop\n");
  for(int iter = 0; iter < niter; ++iter){
    //Rprintf("iteration %i\n", iter);
    // DA step
    //Rprintf("DA step 1");
    for(int k = 0; k < nblocks; ++k){
      rho(k) = R::pnorm(m(k), tau(k), sigma, 0, 0); // = P(t_k < m_k)
      nrej(k) = R::rnbinom(1, rho(k));
    }
    //Rprintf(" 2\n");
    nrejcumsum = cumsum(nrej + 1);
    arma::vec tfull(nrejcumsum(nblocks - 1));
    bsigma = 0;  // running tally, used in sigma and tau steps -- bsigma == oldsqdiff
    for(int k = 0; k < nblocks; ++k){
      tfull(nrejcumsum(k) - nrej(k) - 1) = tobs(k);
      if(nrej(k) > 0){
	// will be slow - calls an R function
	NumericVector tmiss = rtruncnorm(nrej(k), R_NegInf, m(k), tau(k), sigma);
	tfull.subvec(nrejcumsum(k) - nrej(k), nrejcumsum(k) - 1) = as<arma::vec>(tmiss);
      }
      bsigma +=	sum(pow( tfull.subvec(nrejcumsum(k) - nrej(k) - 1, nrejcumsum(k) - 1) - tau(k), 2.0));
    }
    // sigma / sigma^2 step
    //Rprintf("sigma step 1\n");
    omega = arma::randg(1, arma::distr_param(aomega, 1/(bomega + 1/sigma)))[0];
    asigma = (nblocks + sum(nrej) + vsigma)/2;
    sigma2 = 1/arma::randg(1, arma::distr_param(asigma, 1/(bsigma/2 + omega)))[0];
    sigma = sqrt(sigma2);
    // tau step
    //Rprintf("tau step 1");
    zetaprop = zeta + taurwsd*tauchol*arma::randn(nblocks);
    deltaprop = exp(zetaprop);
    tauprop = cumsum(deltaprop);
    tau0prop.subvec(1, nblocks - 1) = tauprop.subvec(0, nblocks - 2);
    psiprop = psifun(tauprop, gam, tau0prop);
    double newsqdiff = 0;
    //Rprintf(" 2");
    for(int k = 0; k < nblocks; ++k){
      newsqdiff += 
	sum(pow(tfull.subvec(nrejcumsum(k) - nrej(k) - 1, nrejcumsum(k) - 1) - tauprop(k), 2.0));
    }
    //Rprintf(" 3");
    lataunum   = -newsqdiff/(2*sigma2) - tauprop(nblocks - 1)/10 + 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop)))) + 
      sum(zetaprop) + sum(psiprop)*log(phi);
    lataudenom = -bsigma/(2*sigma2)     - tau(nblocks - 1)/10 +
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi)))) + 
      sum(zeta) + sum(psi)*log(phi);
    //Rprintf(" 4");
    if(log(arma::randu(1)[0]) < lataunum - lataudenom){
      //      Rcpp::Rcout << "    tau step  " << sum(psi) << "  " << sum(psiprop) << " accept" << std::endl;
      //Rprintf(" 5-0\n");
      delta = deltaprop;
      tau = tauprop;
      psi = psiprop;
      tau0 = tau0prop;
      zeta = zetaprop;
      taurejs(iter,0) = 0;
    } else {
      //      Rcpp::Rcout << "    tau step  " << sum(psi) << "  " << sum(psiprop) << " reject" << std::endl;
      //Rprintf(" 5-1\n");
      taurejs(iter,0) = 1;
    }
    // phi step
    //Rprintf("phi step 1\n");
    phi = R::rbeta(aphi + sum(psi), bphi + N(nblocks));
    // gamma step
    //Rprintf("gamma step 1");
    // gamma_1
    gamprop = gam(0) + gamrwsds(0)*arma::randn(1)[0];
    gamtemp(0) = gamprop;
    gamtemp(1) = gam(1);
    psiprop = psifun(tau, gamtemp, tau0); 
    lanum = sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop)))) + 
      sum(psiprop)*log(phi) - pow(gamprop - mugam(0), 2.0)/(2*sig2gam(0));
    ladenom = sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi)))) + 
      sum(psi)*log(phi) - pow(gam(0) - mugam(0), 2.0)/(2*sig2gam(0));
    if(log(arma::randu(1)[0]) < lanum - ladenom){
      //      Rcpp::Rcout << "    gam1 step " << sum(psi) << "  " << sum(psiprop) << " accept" << std::endl;
      gam(0) = gamprop;
      psi = psiprop;
      gamrejs(iter, 0) = 0;
    } else {
      //      Rcpp::Rcout << "    gam1 step " << sum(psi) << "  " << sum(psiprop) << " reject" << std::endl;
      gamrejs(iter, 0) = 1;
    }
    // gamma_2
    //Rprintf(" 2\n");
    gamprop = gam(1) + gamrwsds(1)*arma::randn(1)[0];
    gamtemp(0) = gam(0);
    gamtemp(1) = gamprop;
    psiprop = psifun(tau, gamtemp, tau0); 
    lanum = sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psiprop)))) + 
      sum(psiprop)*log(phi) - pow(gamprop - mugam(1), 2.0)/(2*sig2gam(1));
    ladenom = sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi + eta)))) - 
      sum(lgamma(Rcpp::as<Rcpp::NumericVector>(wrap(psi)))) + 
      sum(psi)*log(phi) - pow(gam(1) - mugam(1), 2.0)/(2*sig2gam(1));
    if(log(arma::randu(1)[0]) < lanum - ladenom){
      //      Rcpp::Rcout << "    gam2 step " << sum(psi) << "  " << sum(psiprop) << " accept" << std::endl;
      gam(1) = gamprop;
      psi = psiprop;
      gamrejs(iter, 1) = 0;
    } else {
      //      Rcpp::Rcout << "    gam2 step " << sum(psi) << "  " << sum(psiprop) << " reject" << std::endl;
      gamrejs(iter, 1) = 1;
    }
    // d step
    //Rprintf("d step 1\n");
    Nover = N(nblocks) - M(nblocks);
    arma::vec ldvec( (Nover > 0) ? N(nblocks) : M(nblocks) + 1);
    for(int k = 0; k < nblocks; ++k){
      if(nds(k) > 1){
    	if(N(k) == M(k)){
    	  arma::vec lo = rloggamma(nds(k), alpha, 1);
    	  lo = lo - lo.max();
    	  double ldexpsum = log(arma::sum(exp(lo)));
    	  arma::vec ld = log(D(k)) + lo - ldexpsum;
    	  ldvec.subvec(M(k), M(k+1)-1) = ld;
    	} else {
    	  lowlim = (Dbar(k-1) - D(k-1))/D(k);
    	  if(lowlim > 1){
    	    Rprintf("k = %i, lowlim = %f, Dk = %f, Dk-1 = %f\n", k, lowlim, D(k-1), D(k));
    	  }
    	  // calls from R, slow
    	  NumericVector nvtemp = rbetarejC(1, alpha, alpha*(nds(k) - 1), lowlim);  
    	  d1 = nvtemp(0);
    	  arma::vec loother = rloggamma(nds(k)-1, alpha, 1/alpha);
    	  loother = loother - log(arma::sum(exp(loother))) + log(D(k)) + log(1 - d1);
    	  ldvec(M(k))= log(D(k)) + log(d1);
    	  ldvec.subvec(M(k)+1, M(k+1) - 1) = loother;
    	}
      } else {
    	ldvec(M(k)) = log(D(k));
      }
    }
    //Rprintf(" 2");
    if(Nover > 0){
      lowlim = Dbar(nblocks - 1) - D(nblocks - 1);
      prej = R::pgamma(lowlim, alpha, 1/beta, 1, 0);
      arma::vec ldover(Nover);
      if(prej < 0.95){
    	ldover(0) = log(rgammarej(alpha, 1/beta, lowlim)[0]);
      } else {
    	ldover(0) = log(R::qgamma(prej + arma::randu(1)[0]*(1-prej), alpha, 1/beta, 1, 0));
      }
      if(Nover > 1){
    	ldover.subvec(1, Nover-1) = rloggamma(Nover - 1, alpha, 1/beta);
      }
      ldvec.subvec(M(nblocks), N(nblocks) - 1) = ldover;
      mnldvec = mean(ldvec);
      mnexpldvec = mean(exp(ldvec));
      for(int k=0; k < nblocks; ++k){
    	dk1(k) = exp(ldvec(M(k+1)));
      }
    } else {
      ldvec(M(nblocks)) = rloggamma(1, alpha, 1/beta)[0];
      //ldvec(M(nblocks)) = 0;
      mnldvec = mean(ldvec);
      mnexpldvec = mean(exp(ldvec));
      //mnldvec = mean(ldvec.head(ldvec.size()-1));
      // mnexpldvec = mean(exp(ldvec.head(ldvec.size()-1)));
      for(int k=0; k < nblocks; ++k){ //change to nblocks - 1 to go back
    	dk1(k) = exp(ldvec(M(k+1)));
      }
      //dk1(nblocks - 1) = 0;
    }
    //Rprintf(" 3\n");
    // N step
    //Rprintf("N step 1");
    for(int k = 0; k < nblocks; ++k) {
      if(D(k) + dk1(k) <= Dbar(k)){
     	N(k+1) = M(k+1);
      } else {
     	if(M(k+1) >  N(k)){
     	  lower = M(k + 1) - N(k);
     	} else {
     	  lower = 0;
     	}
     	if( k < nblocks - 1){
     	  upper = N(k+2) - N(k);
     	  NumericVector etalogprobs(upper - lower + 1);
     	  NumericVector etaprobs(upper - lower + 1);
     	  NumericVector etaset(upper - lower + 1);
     	  for(int i = 0; i < upper - lower + 1; ++i){
     	    etalogprobs(i) = lgamma(lower + i + psi(k)) - lgamma(lower + i + 1) +
     	      lgamma(upper - lower - i + psi(k+1)) - lgamma(upper - lower - i + 1);
     	    etaset(i) = lower + i;
     	  }
     	  etaprobs = exp(etalogprobs - max(etalogprobs));
     	  etak = sample(etaset, 1, 1, etaprobs)[0];
     	} else {
     	  if(lower > 0){
     	    // calls an R function, will be slow
     	    NumericVector nvtemp = rktnb(1, psi(k), lower - 1, psi(k)*(1-phi)/phi); 
     	    etak = nvtemp(0);
     	  } else {
     	    etak = R::rnbinom(psi(k), phi);
     	  }
     	}
     	N(k + 1) = N(k) + etak;
      }
    }
    //Rprintf(" 2\n");
    eta = diff(N); // update etas after Ns fully sampled
    // Dbar step
    //Rprintf("Dbar step 1\n");
    for(int k = 0; k < nblocks; ++k){
      idxklow = 0;
      while(Dbars(idxklow) < D(k)){
	idxklow += 0;
      }
      idxkhigh = Dbars.size() - 1;
      if(N(k+1) > M(k+1)){
	while(Dbars(idxkhigh) > D(k) + dk1(k)){
	  idxkhigh -= 0;
	}
      }
      arma::vec Dbarproparma = Dbars.subvec(idxklow, idxkhigh);
      NumericVector Dbarprop = Rcpp::as<Rcpp::NumericVector>(wrap(Dbarproparma));
      NumericVector Dbarprob(idxkhigh - idxklow + 1, 1.0/(idxkhigh - idxklow + 1));
      Dbar(k) = sample(Dbarprop, 1, 1, Dbarprob)[0];
    } 
    // alpha step
    //Rprintf("alpha and beta step 1\n");
    alphaprop = exp(log(alpha) + alpharwsd*arma::randn(1)[0]);
      //      double con1 = N(nblocks)*mnldvec - N(nblocks)*log(bbeta + N(nblocks)*mnexpldvec) - balpha;
      //      latau = (alphaprop - alpha)*con1 + (aalpha - 1)*log(alphaprop/alpha) - 
      //      	N(nblocks)*(lgamma(alphaprop) - lgamma(alpha)) +
      //      	lgamma(alphaprop*N(nblocks) + abeta - 1) - lgamma(alpha*N(nblocks) + abeta - 1);
    double con1 = -balpha + N(nblocks)*log(beta) + N(nblocks)*mnldvec;
    double latau = (alphaprop - alpha)*con1 + aalpha*log(alphaprop/alpha) - 
      N(nblocks)*lgamma(alphaprop) + N(nblocks)*lgamma(alpha);
    if(log(arma::randu(1)[0]) < latau){
      alpha = alphaprop;
      alpharejs(iter, 0) = 0;
    } else {
      alpharejs(iter, 0) = 1;
    }
    // beta step
    beta = arma::randg(1, arma::distr_param(abeta + alpha*N(nblocks), scalebeta))[0];   
    // assignment step
    //Rprintf("assignment step 1\n");
    draws(iter,0) = sigma;
    draws(iter,1) = phi;
    draws.submat(iter, 2, iter, 3) = gam.t();
    draws(iter, 4) = alpha;
    draws(iter, 5) = beta;
    draws.submat(iter, 6, iter, nblocks + 5) = 
      arma::conv_to<arma::rowvec>::from(N.subvec(1, nblocks));
    draws.submat(iter, nblocks + 6, iter, 2*nblocks + 5) = tau.t();
    draws.submat(iter, 2*nblocks + 6, iter, 3*nblocks + 5) = Dbar.t();
    // tuning step
    //Rprintf("tuning step 1\n");
    if(tune == 1 && (iter+1) % H == 0){
      Rprintf("iter %i, tuning!\n", iter);
      double acctemp;
      arma::vec meantemp1 = trans(mean(taurejs.submat(iter - H + 1, 0, iter, 0))); 
      acctemp = 1 - meantemp1(0);
      if(acctemp > hightarget(0)){
	taulogrwsd += rwc;
      } else if(acctemp < lowtarget(0)){
	taulogrwsd -= rwc;
      }
      taurwsd = exp(taulogrwsd);
      arma::vec meantemp2 = trans(mean(gamrejs.submat(iter - H + 1, 0, iter, 1))); 
      acctemp = 1 - meantemp2(0);
      if(acctemp > hightarget(1)){
	gamlogrwsds(0) += rwc;
      } else if(acctemp < lowtarget(1)){
	gamlogrwsds(0) -= rwc;
      }
      acctemp = 1 - meantemp2(1);
      if(acctemp > hightarget(1)){
	gamlogrwsds(1) += rwc;
      } else if(acctemp < lowtarget(1)){
	gamlogrwsds(1) -= rwc;
      }
      gamrwsds = exp(gamlogrwsds);
      meantemp1 = mean(alpharejs.submat(iter - H + 1, 0, iter, 0)); 
      acctemp = 1 - meantemp1(0);
      if(acctemp > hightarget(2)){
	alphalogrwsd += rwc;
      } else if(acctemp < lowtarget(2)){
	alphalogrwsd -= rwc;
      }
      alpharwsd = exp(alphalogrwsd);
    }
  }
  Rprintf("Loop finished, collecting 1");
  List initial = List::create(_["sigma"] = sigma, _["phi"] = phi, _["gam"] = gam,
			      _["alpha"] = alpha, _["beta"] = beta, 
			      _["N"] = Rcpp::as<Rcpp::NumericVector>(wrap(N)), 
			      _["tau"] = Rcpp::as<Rcpp::NumericVector>(wrap(tau)), 
			      _["Dbar"] = Rcpp::as<Rcpp::NumericVector>(wrap(Dbar)));
  Rprintf(" 2");
  List rwsds = List::create(_["taulogrwsd"] = taulogrwsd, 
			    _["gamlogrwsds"] = Rcpp::as<Rcpp::NumericVector>(wrap(gamlogrwsds)),
			    _["alphalogrwsd"] = alphalogrwsd);
  Rprintf(" 3");
  List rejs = List::create(_["taurejs"] = Rcpp::as<Rcpp::NumericMatrix>(wrap(taurejs)), 
			   _["gamrejs"] = Rcpp::as<Rcpp::NumericMatrix>(wrap(gamrejs)), 
			   _["alpharejs"] = Rcpp::as<Rcpp::NumericMatrix>(wrap(alpharejs)));
  Rprintf(" 4\n");
  List out = List::create(_["draws"] = draws, _["initial"] = initial,_["rwsds"] = rwsds, 
			  _["rejs"] = rejs);
  Rprintf("collecting finished, outputing\n");
  return out;
}
