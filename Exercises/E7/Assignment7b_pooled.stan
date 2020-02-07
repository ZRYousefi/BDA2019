
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
  real mu_7;
  ypred = normal_rng(mu, sigma);
  mu_7 = normal_rng(mu, sigma);
}

