// this code works for pool, hierarchical, and seperate model. Only the input data should be designed
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of groups
  int<lower=1,upper=K> x[N]; // group indicator
  vector[N] y; //
}
parameters {
  vector[K] theta;        // group means
  real<lower=0> sigma; // individual std
  real<lower=0> tau; // common std
  real mu; //common mean
}
model {
  sigma ~ cauchy(0,4);
  tau ~ cauchy(0,4);
  mu ~ normal(90,15); // weakly informative prior for mu. mean selected close to the mean of data. High variance is assumed. 
  y ~ normal(theta[x], sigma);
  for (i in 1:K){
    theta[i] ~ normal(mu,tau);
  }
}
generated quantities {
 real ypred[K];
 for (i in 1:K){
  ypred[i] = normal_rng (theta[i] , sigma );
 }
}