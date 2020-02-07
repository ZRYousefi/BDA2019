# install.packages("devtools")
# devtools::install_github("avehtari/BDA_course_Aalto",
#                          subdir = "rpackage")
# 
# install.packages("markmyassignment")
# install.packages("aaltobda")

library(ggplot2)
theme_set(theme_minimal())
library(tidyr)

library(markmyassignment)
exercise_path <- "https://github.com/avehtari/BDA_course_Aalto/blob/master/exercises/tests/ex2.yml"
set_assignment(exercise_path)

library(aaltobda)
data("algae")
head(algae)



x <- seq(0,1,length.out=1000)
pp <- dbeta(x, shape1 = 438, shape2 = 544)
plot(x, pp, type='l')
quantile(pp,c(0.025,0.975))
q <- qbeta(pp, shape1 = 438, shape2 = 544)
plot(x,q,type='l')

prior sensitivity page 38

beta_point_est <- function(prior_alpha, prior_beta, data){
  y <- sum(data)
  n <- length(data)
  E_pi = (prior_alpha+y)/(prior_alpha+prior_beta+n)
  return(E_pi)
}

beta_interval <- function(prior_alpha, prior_beta, data, prob){
  y <- sum(data)
  N <- length(data)
  q <- c(qbeta((1-prob)/2, prior_alpha+y, prior_beta+N-y), qbeta(prob+(1-prob)/2, prior_alpha+y, prior_beta+N-y) )
  return(q)
}

beta_low <- function(prior_alpha, prior_beta, data, pi_0){
  y <- sum(data)
  N <- length(data)
  # cumulative distribution function
  cdf <- pbeta(pi_0, prior_alpha+y, prior_beta+N-y)
  return(cdf)
}

algea_test <- c(0, 1, 1, 0, 0, 0)
beta_point_est(prior_alpha = 2, prior_beta = 10, data = algae)
beta_interval(prior_alpha = 2, prior_beta = 10, data = algae, prob = 0.9)
beta_low(prior_alpha = 2, prior_beta = 10, data = algae, pi_0 = 0.2)

mark_my_assignment()


"Uniform prior "
N <- 2
apr <- 0.5
alpha <- N*apr
beta <- N*(1-apr)
post_mean_u = beta_point_est(prior_alpha = alpha , prior_beta = beta, data = algae)
post_int_u = beta_interval(prior_alpha = alpha, prior_beta = beta, data = algae, prob = 0.9)

post_mean_u
post_int_u

"Different beta prior"
N <- c(2, 5, 10, 20, 100, 200, 500)
apr = 0.2
alpha <- N*apr
beta <- N*(1-apr)
post_mean = c()
post_int = matrix(rep(0,2*length(N)), ncol=2)
for(i in 1:length(N)){
  post_mean[i] = beta_point_est(prior_alpha = alpha[i] , prior_beta = beta[i], data = algae)
  post_int[i, ] = beta_interval(prior_alpha = alpha[i], prior_beta = beta[i], data = algae, prob = 0.9)
}

post_mean
post_int


data = algea_test
a <- 2
b <- 10

y <- sum(data)
N <- length(data)
apr <- y/N

# seq creates evenly spaced values
df1 <- data.frame(pi = seq(0.05, 0.8, 0.001)) 

# dbeta computes the posterior density
df1$p <- dbeta(df1$pi, a+y, b+N-y)

# seq creates evenly spaced values from 2.5% quantile
# to 97.5% quantile (i.e., 95% central interval)
# qbeta computes the value for a given quantile given parameters a and b
df2 <- data.frame(pi = seq(qbeta(0, a+y, b+N-y), qbeta(1, a+y, b+N-y), length.out = 1000))
# compute the posterior density
df2$p <- dbeta(df2$pi, a+y, b+N-y)
quantile(df2$pi,c(0.05,0.95))

ggplot(mapping = aes(pi, p)) +
  geom_line(data = df1) +
  # Add a layer of colorized 95% posterior interval
  geom_area(data = df2, aes(fill='1')) +
  # Add the proportion of girl babies in general population
  geom_vline(xintercept = apr, linetype='dotted') +
  # Decorate the plot a little
  labs(title='Prior Beta(2,10) -> Posterior is Beta(438,544)', y = '') +
  scale_y_continuous(expand = c(0, 0.1), breaks = NULL) +
  scale_fill_manual(values = 'lightblue', labels = '90% posterior interval') +
  theme(legend.position = 'bottom', legend.title = element_blank())


df1 <- data.frame(pi = seq(0.05, 0.8, 0.001))
df1$pu <- dbeta(df1$pi, a+1, b+1)


n <- c(2, 20, 200)
apr <- y/N # 2/6

helpref <- function(n, apr, a, b, df)
  cbind(df, pr = dbeta(df$pi, n*apr, n*(1-apr)), po = dbeta(df$pi, n*apr + a, n*(1-apr) + b), n=n)

df2 <- lapply(n, helpref, apr, a, b, df1) %>% do.call(rbind, args = .) %>%
  gather(grp, p, -c(pi, n), factor_key = T)


df2$title <- factor(paste0('alpha/(alpha+beta)=0.160, alpha+beta=',df2$n))
levels(df2$grp) <- c('Posterior with unif prior', 'Prior', 'Posterior')


ggplot(data = df2) +
  geom_line(aes(pi, p, color = grp)) +
  geom_vline(xintercept = 0.16, linetype = 'dotted') +
  facet_wrap(~title, ncol = 1) +
  labs(x = '', y = '') +
  scale_y_continuous(breaks = NULL) +
  theme(legend.position = 'bottom', legend.title = element_blank())





# seq creates evenly spaced values
df1 <- data.frame(theta = seq(0.375, 0.525, 0.001)) 
a <- 438
b <- 544
# dbeta computes the posterior density
df1$p <- dbeta(df1$theta, a, b)


# seq creates evenly spaced values from 2.5% quantile
# to 97.5% quantile (i.e., 95% central interval)
# qbeta computes the value for a given quantile given parameters a and b
df2 <- data.frame(theta = seq(qbeta(0.025, a, b), qbeta(0.975, a, b), length.out = 100))
# compute the posterior density
df2$p <- dbeta(df2$theta, a, b)
quantile(df2$theta,c(0.025,0.975))

ggplot(mapping = aes(theta, p)) +
  geom_line(data = df1) +
  # Add a layer of colorized 95% posterior interval
  geom_area(data = df2, aes(fill='1')) +
  # Add the proportion of girl babies in general population
  geom_vline(xintercept = 0.488, linetype='dotted') +
  # Decorate the plot a little
  labs(title='Uniform prior -> Posterior is Beta(438,544)', y = '') +
  scale_y_continuous(expand = c(0, 0.1), breaks = NULL) +
  scale_fill_manual(values = 'lightblue', labels = '95% posterior interval') +
  theme(legend.position = 'bottom', legend.title = element_blank())