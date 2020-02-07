# install.packages("devtools")
# Sys.setenv(TAR = "/bin/tar")
# devtools::install_github("avehtari/BDA_course_Aalto",
#                          subdir = "rpackage")

library(markmyassignment)
exercise_path <-
  "https://github.com/avehtari/BDA_course_Aalto/blob/master/exercises/tests/ex5.yml"
set_assignment(exercise_path)
# To check your code/functions, just run
# mark_my_assignment()

library(mvtnorm)
library(ggplot2)
theme_set(theme_minimal())
library(aaltobda)
library(gridExtra)

data("bioassay")

data <- bioassay
mu_alpha <- 0
s_alpha <- 2
mu_beta <- 10
s_beta <- 10
rho <- 0.5
s <-  matrix(c(s_alpha^2, rho*s_alpha*s_beta, rho*s_alpha*s_beta, s_beta^2 ), ncol=2)
mu = c(mu_alpha, mu_beta)

"================================   1-a   =================================="
p_log_prior <- function(alpha, beta){
  x = cbind(alpha, beta)
  d <- dmvnorm(x, mean = mu, sigma = s)
  return(log(d))
}

p_log_posterior <- function(alpha, beta, x=bioassay$x, y=bioassay$y, n=bioassay$n){
  p_log_likelihood <- bioassaylp(alpha, beta, x, y, n)
  post <- p_log_prior(alpha, beta) + p_log_likelihood
  return(post)
}

density_ratio <- function(alpha_propose, alpha_previous,
                          beta_propose, beta_previous,
                          x = bioassay$x, y = bioassay$y, n = bioassay$n){
  p1 <- p_log_posterior(alpha_propose, beta_propose, x=bioassay$x, y=bioassay$y, n=bioassay$n)
  p0 <- p_log_posterior(alpha_previous, beta_previous, x=bioassay$x, y=bioassay$y, n=bioassay$n)
  ratio <- exp(p1 - p0) 
  return(ratio)
}
  
density_ratio(alpha_propose = 1.89, alpha_previous = 0.374,
              beta_propose = 24.76, beta_previous = 20.04,
              x = bioassay$x, y = bioassay$y, n = bioassay$n)

density_ratio(alpha_propose = 0.374, alpha_previous = 1.89,
              beta_propose = 20.04, beta_previous = 24.76,
              x = bioassay$x, y = bioassay$y, n = bioassay$n)


"================================  1-b   =================================="

proposal_distribution <- function(param){
  sigma = matrix(c(1, 2, 2, 5), ncol=2)
  return(rnorm(2, mean = param, sd = sigma))
}


metropolis_bioassay <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal <- proposal_distribution(chain[i,])
    r <- density_ratio(alpha_propose = proposal[1], alpha_previous = chain[i,1],
                       beta_propose = proposal[2], beta_previous = chain[i,2],
                       x = bioassay$x, y = bioassay$y, n = bioassay$n)
    if (runif(1) < r){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

# iterations <- 10^4
# startvalue <- c(-5,5)
# chain <- metropolis_bioassay(startvalue, iterations)
# 
# chain <- chain[(1:iterations),]

"================================   2   =================================="

run_in_chains <- function(n_chains, iterations, startvalues){
  m <- n_chains*2
  n <- round(iterations/2/2)
  chains <- array(dim = c(2*n,2, n_chains))
  for(j in 1:n_chains){
    print(paste0('Chain number ',  j))
    chain <- metropolis_bioassay(startvalues[j,], iterations)
    warmUp <- iterations/2
    new_chain <- chain[-(1:warmUp),]
    new_chain <- new_chain[1:warmUp,]
    chains[, ,j] <- new_chain
  }
  return (chains)
}
  
n_chains <- 5
iterations <- 1000000
startvalues <- matrix(c(-1,1,-5,5,-2,2,4,-5,10,-1), byrow = T, ncol=2)
chains <- run_in_chains(n_chains, iterations, startvalues)
warmUp <- iterations/2
# new_chain <- chain[-(1:warmUp),]

inds <- 1:warmUp
dfs <- data.frame(iter = inds, chains[inds, 1, ]) %>%
  gather(chain, theta1, -iter) %>%
  within(theta2 <- c(chains[inds, 2, ])) %>%
  gather(var, val, -iter, -chain)


#' Plot trends of the all draws
ggplot(data = dfs) +
  geom_line(aes(iter, val, color = chain)) +
  facet_grid(var~.) +
  labs(title = 'Visually converged', y = '') +
  scale_color_discrete(guide = FALSE)

"================================   3   =================================="

R_hat <- function(chains){
  m <- 2*dim(chains)[3]
  n <- round(dim(chains)[1]/2)
  psi_ij <- array(dim = c(n, 2, m))
  for(j in 1:n_chains){    
    psi_ij[, , 2*j-1] <- chains[(1:n), , j]
    psi_ij[, , 2*j] <- chains[-(1:n), , j]
  }
  
  psi_j <- array(dim = c(m, 2))
  s_j_2 <- array(dim = c(m, 2))
  for(j in 1:m){
    psi_j[j,] <- colMeans(psi_ij[, ,j])
    s_j_2[j,] <- rowSums(apply(psi_ij[, ,j], 1, function(x) x-psi_j[j,])^2) / (n-1)
  }
  psi <- colMeans(psi_j)
  
  B <- ( n/(m-1) ) * sum(apply(psi_j, 1, function(x) x-psi)^2)
  W <- colMeans(s_j_2)
  var_hat <- ((n-1)/n) * W  + (1/n) * B 
  R_hat <- sqrt(var_hat/W)
  return(R_hat)
}
  
R_hat(chains)
"================================   4   =================================="


samp_A <- chains[, 1, 1]
samp_B <- chains[, 2, 1]

xl <- c(-5, 10)
yl <- c(-10, 40)

ggplot(data = data.frame(samp_A, samp_B)) +
  geom_point(aes(samp_A, samp_B), color = 'blue', size = 0.3) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta')


