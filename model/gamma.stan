data {
  int<lower=0> nobs;
  vector<lower=0>[nobs] mb;
  vector<lower=0>[2] logalphapars;
  vector<lower=0>[2] betapars;
}
parameters {
  real logalpha;
  real<lower=0> beta;
}
transformed parameters{
  real<lower=0> alpha;
  alpha <- exp(logalpha);
}
model {
  mb ~ gamma(alpha, beta);
  logalpha ~ normal(logalphapars[1], logalphapars[2]); // implies a lognormal prior on alpha
  beta ~ gamma(betapars[1], betapars[2]);
}
