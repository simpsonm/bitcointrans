data {
  int<lower=0> nobs;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
  vector<lower=0>[nobs] mb;
  vector<lower=0>[2] sigpars;
  vector<lower=0>[2] sigepspars;
  vector<lower=0>[2] bpars;
  real<lower=0> omega;
  vector[3] alphamns;
  vector[3] alphasds;
  vector<lower=0>[2] rhopars;
  real<lower=0> phipar;
  vector[2] lgam0pars;
}
parameters {
  vector<lower=0>[nobs] delta;
  real<lower=0> b;
  real<lower=0> sigma;
  vector[nobs+1] lgam_raw;
  real<lower=0> phi;
  real<lower=0> sigeps;
  real alpha0;
  real alpha1;
  real alpha2;
  real<lower=0, upper=1> rho;
}
transformed parameters{
  vector[nobs] tau;
  vector[nobs+1] lgam;
  vector[nobs+1] gam;
  vector[nobs+1]  mu;
  vector[nobs] lgammns;
  vector<lower=0>[nobs] lgamsds;
  real<lower=0> psi;

  tau <- cumulative_sum(delta);

  psi <- exp(-rho);
  mu[1] <- alpha0 + alpha2*sin(phi);
  lgam[1] <- lgam0pars[1] + lgam0pars[2]*lgam_raw[1];
  for(i in 1:nobs){
    mu[i+1] <- alpha0 + alpha1*tau[i] + alpha2*sin(omega*tau[i] + phi);
    lgammns[i] <- exp(-psi*delta[i])*(lgam[i] - mu[i]) + mu[i+1];
    lgamsds[i] <- sigeps*sqrt(1 - exp(-2*psi*delta[i]))/sqrt(2*psi);
    lgam[i+1] <- lgammns[i] + lgam_raw[i+1]*lgamsds[i];
  }
  gam <- exp(lgam);
}
model {
  for(i in 1:nobs){
    delta[i] ~ exponential(lambda);
    t[i] ~ lognormal(log(tau[i]) - sigma^2/2, sigma) T[fmax(lb[i],0),]; 
    mb[i] ~ gamma(b*(gam[i+1] - gam[i]), b);
    lgam_raw[i] ~ normal(0, 1);
  }
  sigma ~ normal(log(sigpars[1]), sigpars[2]);
  b ~ gamma(bpars[1], bpars[2]);
  sigeps ~ cauchy(log(sigepspars[1]), sigepspars[2]);
  alpha0 ~ normal(alphamns[1], alphasds[1]);
  alpha1 ~ normal(alphamns[2], alphasds[2]);
  alpha2 ~ normal(alphamns[3], alphasds[3]);
  rho ~ beta(rhopars[1], rhopars[2]);
  phi ~ uniform(0, phipar);
}
