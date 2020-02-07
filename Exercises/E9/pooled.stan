
data {
  int<lower=0> N; // number of data points
  vector[N] y; //
}
parameters {
  real mu; // group means
  real<lower=0> sigma; // common std
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
  real ypred;
  vector[N] log_lik;
  ypred = normal_rng(mu, sigma);
  for (i in 1:N){
    log_lik[i] = normal_lpdf(y[i] | mu, sigma);
  }
}

