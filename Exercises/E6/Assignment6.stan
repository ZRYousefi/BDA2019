// bioassay model
data {
  int<lower=0> k;
  int n[k];
  row_vector[k] x;  
  int y[k]; 
  vector[2] mu_0;
  matrix[2, 2] sigma_0;
}

parameters {
  row_vector[2] theta;
}

model {
  theta ~ multi_normal(mu_0, sigma_0);
  for (i in 1:k) {
    y[i] ~ binomial_logit(n[i],theta[1]+theta[2]*x[i]);
  }
}


