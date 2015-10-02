data {
  int<lower=0> nobs;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
  vector<lower=0>[nobs] mb;
  vector<lower=0>[2] sigpars;
  vector<lower=0>[2] apars;
  vector<lower=0>[2] bpars;
}
parameters {
  vector<lower=0>[nobs] x;
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> sigma;
}
transformed parameters{
  vector[nobs] tau;
  real<lower=0> mu;
  tau <- cumulative_sum(x);
  mu <- a/b;
}
model {
  for(i in 1:nobs){
    x[i] ~ exponential(lambda);
    t[i] ~ lognormal(log(tau[i]) - sigma^2/2, sigma) T[fmax(lb[i],0),]; 
    mb[i] ~ gamma(a*x[i], b);
  }
  sigma ~ cauchy(log(sigpars[1]), sigpars[2]);
  a ~ gamma(apars[1], apars[2]);
  b ~ gamma(bpars[1], bpars[2]);
}
