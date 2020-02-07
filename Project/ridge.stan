data{
  int N_train;             // training observations
  int N_test;              // test observations
  int N_pred;              // predictor variables
  vector[N_train] y_train; // training outcomes
  matrix[N_train, N_pred] X_train; // training data
  matrix[N_test, N_pred] X_test;   // testing data
}
parameters{
  real beta_0; // intercept
  real<lower=0> sigma;   // error SD
  real<lower=0> sigma_B; // SD across betas
  vector[N_pred] beta;   // regression beta weights
}
model{
  beta_0 ~ normal(0, 1); // prior on intercept
  sigma ~ normal(0, 1); // model error SD
  beta ~ normal(0, sigma_B); // beta prior : ridge 
  y_train ~ normal(beta_0 + X_train * beta, sigma); // model likelihood
}
generated quantities{ 
  real y_test[N_test]; // test data predictions
  vector[N_train] log_lik;
  for(i in 1:N_train){
  log_lik[i] = normal_lpdf(y_train[i] | beta_0 + X_train[i,] * beta, sigma);
  }
  for(i in 1:N_test){
      y_test[i] = normal_rng(beta_0 + X_test[i,] * beta, sigma);
  }
}


