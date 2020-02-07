data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  real mu0; // prior mean
  real<lower=0> sigma0; // prior std
  vector[K] mu; // group means
  real<lower=0> sigma; // group stds
}
model {
  mu ~ normal(mu0, sigma0); // population prior with unknown parameters
  y ~ normal(mu[x] , sigma);
}
generated quantities {
  real ypred[K];
  for (i in 1:K){
    ypred[i] = normal_rng(mu[i], sigma);
  }
}

