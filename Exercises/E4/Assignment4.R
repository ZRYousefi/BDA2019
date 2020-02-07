# install.packages("devtools")
# Sys.setenv(TAR = "/bin/tar")
# devtools::install_github("avehtari/BDA_course_Aalto",
#                          subdir = "rpackage")

library(markmyassignment)
exercise_path <- "https://github.com/avehtari/BDA_course_Aalto/blob/master/exercises/tests/ex4.yml"
set_assignment(exercise_path)
# To check your code/functions, just run
mark_my_assignment()

#install.packages("mvtnorm")

library(ggplot2)
theme_set(theme_minimal())
library(aaltobda)


data("bioassay")

data <- bioassay
mu_alpha <- 0
s_alpha <- 2
mu_beta <- 10
s_beta <- 10
rho <- 0.5


p_log_prior <- function(alpha, beta){
  x = cbind(alpha, beta)
  d <- dmvnorm(x, mean=c(mu_alpha, mu_beta), 
               sigma = matrix(c(s_alpha^2, rho*s_alpha*s_beta, rho*s_alpha*s_beta, s_beta^2 ), ncol=2))
  return(log(d))
}

# logit_inv <- function(theta){
#   return (1/(1+exp(-theta)))
# }


p_log_posterior <- function(alpha, beta, x=bioassay$x, y=bioassay$y, n=bioassay$n){
  # ll <- c()
  # p_log_likelihood <- 0
  # for(i in 1:length(x)){
  #   ll[i] <- y[i] * log(logit_inv(alpha+beta*x[i])) + (n[i]-y[i]) * log(1-logit_inv(alpha+beta*x[i]))
  #   p_log_likelihood <- p_log_likelihood + ll[i]
  #   # ll[i] <- (logit_inv(alpha+beta*x[i]))^y[i] * (1-logit_inv(alpha+beta*x[i]))^(n[i]-y[i])
  # }
  # p_log_likelihood <- log(Reduce(`*`,ll))
  p_log_likelihood <- bioassaylp(alpha, beta, x, y, n)
  post <- p_log_prior(alpha, beta) + p_log_likelihood
  return(post)
}


alpha <- 3
beta <- 9
p_log_prior(alpha,beta)
p_log_posterior(alpha,beta, x=bioassay$x, y=bioassay$y, n=bioassay$n)

mark_my_assignment()

bioassay_posterior_density_plot(alpha_limits = c(-4, 4),
                                beta_limits = c(-10, 30), x=bioassay$x, y=bioassay$y, n=bioassay$n) 


"================  e    ================"


log_importance_weights <- function(alpha, beta){
  S <- length(alpha)
  w <- c()
  for(i in 1:S){
    w[i] <-  p_log_posterior(alpha[i], beta[i], x=bioassay$x, y=bioassay$y, n=bioassay$n) - p_log_prior(alpha[i], beta[i]) 
}  
  return(w)
  # return(round(w, digits = 2))
} 

normalized_importance_weights <- function(alpha, beta){
  log_w <- log_importance_weights(alpha, beta)
  exp_w <- exp(log_w)
  return(exp_w/sum(exp_w))
  #return(round(exp_w/sum(exp_w), digits = 3))
}

alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
log_importance_weights(alpha, beta)
# [1] -8.95 -23.47 -6.02 -8.13 -16.61 -14.57
normalized_importance_weights(alpha, beta)
# [1] 0.045 0.000 0.852 0.103 0.000 0.000

mark_my_assignment()

"===============   f   ================="


nr <-5000
r <- rmvnorm(nr, mean=c(mu_alpha, mu_beta),
             sigma = matrix(c(s_alpha^2, rho*s_alpha*s_beta,
                              rho*s_alpha*s_beta, s_beta^2 ), ncol=2))
alpha <- r[, 1]
beta <- r[, 2]

posterior_mean(alpha, beta)
S_eff(alpha, beta)
# alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
# beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
# 
# alpha = c(1.4, 1.5, -2.3, -2.3, -0.5, 1.9, 1.6, -1.8, -2.7, -1.9)
# beta = c(11.4, 8.2, 10.8, 7.4, 18.2, 27.8, 5.1, 1.9, -15.7, 25.1)

posterior_mean <- function(alpha, beta){
  ab <- cbind(alpha, beta)
  colnames(ab) <- NULL
  post <- colSums( ab * normalized_importance_weights(alpha, beta) )
  #return(round(post, digit=3))
  return(post)
}

posterior_mean(alpha, beta)

"==================   g   ========================"
S_eff <- function(alpha, beta){
  s_eff <- 1/sum(normalized_importance_weights(alpha, beta)^2)
  return (s_eff)
}

S_eff(alpha, beta)

mark_my_assignment()


"==================   h   ========================="
# A = seq(-2, 6, length.out = 100)
# B = seq(-2, 30, length.out = 100)
# make vectors that contain all pairwise combinations of A and B
# cA <- rep(A, each = length(B))
# cB <- rep(B, length(A))


nr <- 5000
cR <- rmvnorm(nr, mean=c(mu_alpha, mu_beta),
             sigma = matrix(c(s_alpha^2, rho*s_alpha*s_beta,
                              rho*s_alpha*s_beta, s_beta^2 ), ncol=2)) 
cA <- cR[, 1]
cB <- cR[, 2]


nsamp <- 1000
samp_indices <- sample(length(cA), size = nsamp, replace = FALSE, prob = normalized_importance_weights(cA, cB) )


samp_A <- cA[samp_indices[1:nsamp]]
samp_B <- cB[samp_indices[1:nsamp]]
# add random jitter, see BDA3 p. 76
# samp_A <- samp_A + runif(nsamp, A[1] - A[2], A[2] - A[1])
# samp_B <- samp_B + runif(nsamp, B[1] - B[2], B[2] - B[1])


xl <- c(-5, 5)
yl <- c(-11, 31)

ggplot(data = data.frame(samp_A, samp_B)) +
  geom_point(aes(samp_A, samp_B), color = 'blue', size = 0.3) +
  coord_cartesian(xlim = xl, ylim = yl) +
  labs(x = 'alpha', y = 'beta')

"===================   i   ===================="
bpi <- samp_B > 0
samp_ld50 <- -samp_A[bpi]/samp_B[bpi]

p_positive_beta <- length(bpi)/nsamp

"==================  j   ======================"
his <- ggplot() +
  geom_histogram(aes(samp_ld50), binwidth = 0.04,
                 fill = 'steelblue', color = 'black') +
  coord_cartesian(xlim = c(-0.8, 0.8)) +
  labs(x = 'LD50 = -alpha/beta')

grid.arrange(sam, his, ncol = 2)
