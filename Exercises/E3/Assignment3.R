install.packages("devtools")
devtools::install_github("avehtari/BDA_course_Aalto",
                         subdir = "rpackage")

install.packages("markmyassignment")
install.packages("aaltobda")

library(ggplot2)
theme_set(theme_minimal())
library(tidyr)

library(markmyassignment)
exercise_path <- "https://github.com/avehtari/BDA_course_Aalto/blob/master/exercises/tests/ex3.yml"
set_assignment(exercise_path)


library(aaltobda)


"====== noninformative prior ========"

data("windshieldy1")
head(windshieldy1)

data <- windshieldy1
n <- length(data)
mu <- mean(data)
sigma <- sd(data) 


x <- seq(12, 18, length.out = 1000)
exact_posterior_mu <- dtnew(x, df=n-1, mean=mu, scale=sigma/sqrt(n)) # dt((x-mu) / sigma/sqrt(n), df=n-1) / sigma/sqrt(n)
emprical_posterior_mu <- dnorm(x, mu, sigma/sqrt(n))


mu_point_est <- function(data){
  n <- length(data)
  mu <- mean(data)
  sigma <- sd(data)
  rtg <- rt(100000, df=n-1)
  rr <- (rtg * sigma/sqrt(n) ) + mu
  mu_post <- mean(rr)

  return(mu_post)
}  

mu_interval <- function(data, prob){
  n <- length(data)
  mu <- mean(data)
  sigma <- sd(data)
  rtg <- rt(100000, df=n-1)
  rr <- (rtg * sigma/sqrt(n) ) + mu
  q <- quantile(rr, c((1-prob)/2, prob+(1-prob)/2), names = FALSE)
  return(q)
}

mu_point_est(data)
mu_interval(data, prob = 0.95)

ggplot() +
  geom_line(aes(x, exact_posterior_mu, color='exact')) +
  geom_line(aes(x, emprical_posterior_mu, color='emprical')) + 
  geom_vline(aes(xintercept = mu_point_est(data), color = 'posterior mean'),
             linetype = 'dashed', show.legend = F) +
  geom_vline(aes(xintercept = c(mu_interval(data, prob = 0.95)), color = '95% interval'),
             linetype = 'solid', show.legend = F)
  labs(title = 'Marginal of mu', x = 'mu', y = '')
  


windshieldy_test <- c(13.357, 14.928, 14.896, 14.820)





"====== predictive distribution ==========="



x <- seq(5, 20, length.out = 1000)
exact_posterior_pred <- dtnew(x, df=n-1, mean=mu, scale=sqrt(sigma*sqrt(1+1/n))) # dt((x-mu) / sigma*sqrt(1+1/n), n-1) / sigma*sqrt(1+1/n)
emprical_posterior_pred <- dnorm(x, mu, sqrt(sigma*sqrt(1+1/n)))


mu_pred_point_est <- function(data){
  n <- length(data)
  mu <- mean(data)
  sigma <- sd(data)
  rtg <- rt(100000, df=n-1)
  rr <- (rtg * sigma*sqrt(1+(1/n)) ) + mu
  mu_post <- mean(rr)
  return(mu_post)
}

mu_pred_interval <- function(data, prob = 0.95){
  n <- length(data)
  mu <- mean(data)
  sigma <- sd(data)
  rtg <- rt(100000, df=n-1)
  rr <- (rtg * sigma*sqrt(1+(1/n)) ) + mu
  q <- quantile(rr, c((1-prob)/2, prob+(1-prob)/2))
  return(q)
}

mu_pred_point_est(data)
mu_pred_interval(data, prob = 0.95)




ggplot() +
  geom_line(aes(x, exact_posterior_pred, color='exact')) +
  geom_line(aes(x, emprical_posterior_pred, color='emprical')) + 
  geom_vline(aes(xintercept = mu_pred_point_est(data), color = 'posterior mean'),
             linetype = 'dashed', show.legend = F) +
  geom_vline(aes(xintercept = c(mu_pred_interval(data, prob = 0.95)), color = '95% interval'),
             linetype = 'solid', show.legend = F) +
labs(title = 'Marginal of predicted mu', x = 'mu', y = '')






"==================  2 ==================="


n0 <- 674
y0 <- 39
n1 <- 680
y1 <- 22
a0 <- 1
b0 <- 1
a1 <- 1
b1 <- 1

post_alpha0 <- a0 + y0
post_beta0 <- b0 + n0 - y0
post_dist0 <- rbeta(1000, post_alpha0, post_beta0)

post_alpha1 <- a1 + y1
post_beta1 <- b1 + n1 - y1
post_dist1 <- rbeta(1000, post_alpha1, post_beta1)



posterior_odds_ratio_point_est <- function(p0, p1){
  psi <- (p1/(1-p1))/(p0/(1-p0))
  return(mean(psi))
}

posterior_odds_ratio_interval <- function(p0, p1, prob = 0.9){
  psi <- (p1/(1-p1))/(p0/(1-p0))
  q <- c(quantile(psi, (1-prob)/2), quantile(psi, prob+(1-prob)/2))
  return(q)
}

odds_ratio <- (post_dist1/(1-post_dist1))/(post_dist0/(1-post_dist0))

ggplot() +
  geom_histogram(aes(odds_ratio), binwidth = 0.09, fill = 'steelblue', color = 'black') +
  coord_cartesian(xlim = c(0, 1.5)) +
  scale_y_continuous(breaks = NULL) +
  labs(title = 'Odds ratio histogram', x = 'odss_ratio')+
  geom_vline(aes(xintercept = mean(odds_ratio), color = 'q'),
             linetype = 'dashed', show.legend = F)

posterior_odds_ratio_point_est(post_dist0, post_dist1)
posterior_odds_ratio_interval(post_dist0, post_dist1, prob = 0.95)


A0 <- c(1, 2, 0.5, 5)
B0 <- c(1, 10, 10 ,100)
A1 <- c(1, 2, 0.4, 4)
B1 <- c(1, 10, 10, 100)
post_mean = c()
post_int = matrix(rep(0,2*length(A0)), ncol=2)
for(i in 1:length(A0)){
  a0 <- A0[i]
  b0 <- B0[i]
  a1 <- A1[i]
  b1 <- B1[i]
  
  post_alpha0 <- a0 + y0
  post_beta0 <- b0 + n0 - y0
  prior_dist0 <- rbeta(1000, a0, b0)
  post_dist0 <- rbeta(1000, post_alpha0, post_beta0)
  
  post_alpha1 <- a1 + y1
  post_beta1 <- b1 + n1 - y1
  prior_dist1 <- rbeta(1000, a1, b1)
  post_dist1 <- rbeta(1000, post_alpha1, post_beta1)
  
  prior_dist <- rbeta(1000, a0+a1, b0+b1)
  post_mean[i] <- posterior_odds_ratio_point_est(post_dist0, post_dist1)
  post_int[i, ] <- posterior_odds_ratio_interval(post_dist0, post_dist1, prob = 0.95)
  
  
}

post_mean
post_int

set.seed(4711)
p0 <- rbeta(100000, 5, 95)
p1 <- rbeta(100000, 10, 90)
posterior_odds_ratio_point_est(p0, p1)
posterior_odds_ratio_interval(p0, p1, prob = 0.9)

mark_my_assignment()



"========================   3   ==========================="

data("windshieldy1")
data("windshieldy2")


post_mean <- function(data){
  n <- length(data)
  mu <- mean(data)
  sigma <- sd(data)
  rtg <- rt(100000, df=n-1)
  rr <- (rtg * sigma/sqrt(n) ) + mu
  return(rr)
}


  
data1 <- windshieldy1
data2 <- windshieldy2 

mu_difference <- post_mean(data1) - post_mean(data2)
posterior_mean <- mean(mu_difference)
posterior_mean

prob <- 0.95
posterior_interval<- quantile(mu_difference, c((1-prob)/2, prob+(1-prob)/2), names = FALSE)
posterior_interval


labs <- c('posterior mean')
ggplot() +
  geom_histogram(aes(mu_difference), binwidth = 0.6, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of mu_d ', x = 'mu_d')+
  geom_vline(aes(xintercept = posterior_mean, color = 'q'),
             linetype = 'dashed', show.legend = F) + 
  geom_vline(aes(xintercept = c(posterior_interval), color = '95% interval'),
             linetype = 'solid', show.legend = F) 

mark_my_assignment()
