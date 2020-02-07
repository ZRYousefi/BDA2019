data{
    int N_train;             // "# training observations"
    int N_test;              // "# test observations"
    int N_pred;              // "# predictor variables"
    vector[N_train] y_train; // "training outcomes"
    matrix[N_train, N_pred] X_train; // "training data"
    matrix[N_test, N_pred] X_test;   // "testing data"
}
parameters{
	real beta_0; //intercept
	real<lower=0> sigma2; //error variance
	vector[N_pred] beta_raw; // regression parameters
	//hyperparameters prior
	vector<lower=0, upper=pi()/2>[N_pred] tau_unif;
	real<lower=0> lambda;
}
transformed parameters{
	vector[N_pred] beta;
	real<lower=0> sigma; //error sd
	vector<lower=0>[N_pred] tau;
	vector[N_train] linpred; //mean normal model
	tau = lambda * tan(tau_unif); // implies tau ~ cauchy(0, lambda)
	for(j in 1:N_pred){
		beta[j] = tau[j] * beta_raw[j];
	}
	sigma = sqrt(sigma2);
	linpred = beta_0 + X_train*beta;
}
model{
 //prior regression coefficients: horseshoe
	beta_raw ~ normal(0, 1); //implies beta ~ normal(0, tau)
	
	lambda ~ cauchy(0, 10000);
	
 //priors nuisance parameters: uniform on log(sigma^2) & beta_0
	target += -2 * log(sigma); 
	
 //likelihood
	y_train ~ normal(linpred, sigma);
}
generated quantities{ //predict responses test set
	real y_test[N_test]; //predicted responses
	vector[N_train] log_lik;
	for(i in 1:N_train){
    log_lik[i] = normal_lpdf(y_train[i] | beta_0 + X_train[i,] * beta, sigma);
  }
	for(i in 1:N_test){
		y_test[i] = normal_rng(beta_0 + X_test[i,] * beta, sigma);
	}
}	
	
	
	