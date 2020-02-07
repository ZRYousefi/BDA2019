data {
  int<lower=0> N; // number of data points
  vector[N] x; // observation year
  vector[N] y; // observation number of drowned
  real xpred; // prediction year
 }
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
transformed parameters {
 vector[N] mu;
   mu = alpha + beta *x;
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
  real ypred;
  ypred = normal_rng(alpha + beta*xpred, sigma);
}



