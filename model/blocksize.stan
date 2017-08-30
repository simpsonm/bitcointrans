data {
  int<lower=1> n_blocks;
  int<lower = 1> n_betas;
  vector<lower = 0, upper = 2>[n_betas] alpha; // polynomial degrees in regression
  vector[n_blocks] t;                          // arrival time for each block
  vector[n_blocks] lb;                         // lower bound for arrival time of each block
  real<lower=0> lambda;                        // 1 / mean arrival time, typically 1/10
  vector<lower=0>[n_blocks] d;                 // block sizes (in MBs)
  // prior hyperparameters
  real<lower = 0> sigma_prior_scale;
  real beta_prior_mean;
  real<lower = 0> beta_prior_scale;
  real<lower = 0> theta_prior_scale;
}
parameters {
  vector<lower=0>[n_blocks] delta;        // inter-arrival times
  real<lower=0> sigma;                    // SD of arrival time measurement error
  vector<lower = 0>[n_betas] beta;        // vector of regression coefficients
  real<lower = 0> theta;                  // block size dispersion parameter
 
}
transformed parameters{
  vector[n_blocks] tau;                   // arrival times
  vector[n_blocks] psi;                   // integrated intensities
  vector[n_beta] beta_over_alpha;         // beta_p/alpha_p for all p

  tau = cumulative_sum(delta);
  beta_over_alpha = beta ./ alpha;
  for(i in 1:n_blocks){
    psi[i] = dot_product(beta_over_alpha, tau[i]^alpha);
    if(i > 1){
      psi[i] = psi[i] - psi[i-1];   // psi_k = sum_p beta_p / alpha_p * (tau_{k}^{alpha_p} - tau_{k-1}^{alpha_p})
    }
  }
}
model {
  delta ~ exponential(lambda);                      // inter-arrival time model
  d ~ gamma(theta*psi, theta);                      // block size model
  for(i in 1:n_blocks){
    t[i] ~ normal(tau[i], sigma) T[fmax(lb[i],0),]; // observed time measurement error model
  }
  // priors
  sigma ~ cauchy(0, sigma_prior_scale);
  beta ~ normal(beta_prior_mean, beta_prior_scale);
  theta ~ cauchy(0, theta_prior_scale);
}
