data {
  int<lower=0> nobs;
  //  int<lower=0> nprev;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
}
parameters {
  vector<lower=0.01>[nobs] x;
  real<lower=0.01> sigma;
}
transformed parameters{
  //  positive_ordered[nobs] tau;
  vector[nobs] tau;
  //  vector[nobs] alpha;
  //  vector[nobs] beta;
  tau <- cumulative_sum(x);
  //  beta <- tau / sigma^2;
  //  alpha <- beta .* tau;
}
model {
  //  t[1] ~ normal(tau[1], sigma) T[lb[1],]; 
  //  x[1] ~ gamma(nprev, lambda);  // uses shape/rate parameterization
  for(i in 1:nobs){
    t[i] ~ lognormal(log(tau[i] - sigma^2/2), sigma) T[fmax(lb[i],0),]; 
    x[i] ~ exponential(lambda);
  }
  sigma ~ cauchy(0, 2.5);
}
