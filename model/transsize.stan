data {
  int<lower=1> n_blocks;
  vector<lower=0>[n_blocks] d;                 // block sizes (in MBs)
  int<lower=0>[n_blocks] n_trans;              // transaction counts
  // priors
  real<lower = 0> alpha_prior_scale;
  real<lower = 0> beta_prior_scale;
}
parameters {
  real<lower = 0> alpha;
  real<lower = 0> beta;
}
model {
  for(i in 1:n_blocks){
    d[i] ~ gamma(alpha*n_trans[i], beta);
  }
  alpha ~ cauchy(0, alpha_prior_scale);
  beta ~ cauchy(0, beta_prior_scale);
}
