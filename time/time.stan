data {
  int<lower=0> nobs;
  int<lower=0> nprev;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
}
parameters {
  vector<lower=0>[nobs] x;
  real<lower=0> sigma;
}
transformed parameters{
  //  positive_ordered[nobs] tau;
  vector[nobs] tau;
  tau <- cumulative_sum(x);
}
model {
  t[1] ~ normal(tau[1], sigma) T[lb[1],]; 
  x[1] ~ gamma(nprev, lambda);  // uses shape/rate parameterization
  for(i in 2:nobs){
    t[i] ~ normal(tau[i], sigma) T[lb[i],]; 
    x[i] ~ exponential(lambda);
  }
  sigma ~ cauchy(0, 2.5) T[0,];
}
