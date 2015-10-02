// doesn't work due to numerical issues

data {
  int<lower=0> nobs;
  vector[nobs] t;
  vector[nobs] lb;
  real<lower=0> lambda;  // 1 / mean arrival time
  vector<lower=0>[nobs] mb;
  vector<lower=0>[2] betapars;
  vector<lower=0>[2] apars;
  vector<lower=0>[2] bpars;

}
parameters {
  vector<lower=0>[nobs] x;
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> beta;
}
transformed parameters{
  vector[nobs] tau;
  real<lower=0> mu;
  vector<lower=0>[nobs] alpha;
  tau <- cumulative_sum(x);
  mu <- a/b;
  alpha <- beta*tau;
}
model {
  for(i in 1:nobs){
    x[i] ~ exponential(lambda);
    t[i] ~ gamma(alpha[i], beta) T[fmax(lb[i],0),]; 
    mb[i] ~ gamma(a*x[i], b);
  }
  beta ~ gamma(betapars[1], betapars[2]);
  a ~ gamma(apars[1], apars[2]);
  b ~ gamma(bpars[1], bpars[2]);
}
