data {
  int<lower=0> nobs;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
  vector<lower=0>[nobs] x;
  vector[2] sigpars;
  vector<lower=0>[2] loggammapars;
  vector<lower=0>[2] betapars;
}
parameters {
  vector<lower=0>[nobs] delta;
  real loggamma;
  real<lower=0> beta;
  real<lower=0> sigma;
}
transformed parameters{
  vector[nobs] tau;
  real<lower=0> gamma;
  tau <- cumulative_sum(delta);
  gamma <- exp(loggamma);
}
model {
  for(i in 1:nobs){
    delta[i] ~ exponential(lambda);
    t[i] ~ lognormal(log(tau[i]) - sigma^2/2, sigma) T[fmax(lb[i],0),]; 
    x[i] ~ gamma(beta * gamma * delta[i], beta);
  }
  sigma ~ cauchy(sigpars[1], sigpars[2]);
  loggamma ~ normal(loggammapars[1], loggammapars[2]);
  beta ~ gamma(betapars[1], betapars[2]);
}
