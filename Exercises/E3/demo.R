library(ggplot2)
theme_set(theme_minimal())
library(grid)
library(gridExtra)
library(tidyr)



"Data"

y <- c(93, 112, 122, 135, 122, 150, 118, 90, 124, 114)

"Sufficient statistics"

n <- length(y)
s2 <- var(y)
my <- mean(y)

"Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y) Sample from the joint posterior using this factorization"

# helper functions to sample from and evaluate
# scaled inverse chi-squared distribution
rsinvchisq <- function(n, nu, s2, ...) nu*s2 / rchisq(n , nu, ...)
dsinvchisq <- function(x, nu, s2){
  exp(log(nu/2)*nu/2 - lgamma(nu/2) + log(s2)/2*nu - log(x)*(nu/2+1) - (nu*s2/2)/x)
}

"Sample 1000 random numbers from p(sigma2|y)"

ns <- 1000
sigma2  <- rsinvchisq(ns, n-1, s2)

"Sample from p(mu|sigma2,y)"

mu <- my + sqrt(sigma2/n)*rnorm(length(sigma2))

"Create a variable sigma and sample from predictive distribution p(ynew|y) for each draw of (mu, sigma)"

sigma <- sqrt(sigma2)
ynew <- rnorm(ns, mu, sigma)

"For mu, sigma and ynew compute the density in a grid ranges for the grids"

t1l <- c(90, 150)
t2l <- c(10, 60)
nl <- c(50, 185)
t1 <- seq(t1l[1], t1l[2], length.out = ns)
t2 <- seq(t2l[1], t2l[2], length.out = ns)
xynew <- seq(nl[1], nl[2], length.out = ns)

"Compute the exact marginal density of mu"

# multiplication by 1./sqrt(s2/n) is due to the transformation of
# variable z=(x-mean(y))/sqrt(s2/n), see BDA3 p. 21
pm <- dt((t1-my) / sqrt(s2/n), n-1) / sqrt(s2/n)

"Estimate the marginal density using samples and ad hoc Gaussian kernel approximation"

pmk <- density(mu, adjust = 2, n = ns, from = t1l[1], to = t1l[2])$y

"Compute the exact marginal density of sigma"

# the multiplication by 2*t2 is due to the transformation of
# variable z=t2^2, see BDA3 p. 21
ps <- dsinvchisq(t2^2, n-1, s2) * 2*t2

"Estimate the marginal density using samples and ad hoc Gaussian kernel approximation"

psk <- density(sigma, n = ns, from = t2l[1], to = t2l[2])$y

"Compute the exact predictive density"

# multiplication by 1./sqrt(s2/n) is due to the transformation of variable
# see BDA3 p. 21
p_new <- dt((xynew-my) / sqrt(s2*(1+1/n)), n-1) / sqrt(s2*(1+1/n))

"Evaluate the joint density in a grid. Note that the following is not normalized, but for plotting contours it does not matter."

# Combine grid points into another data frame
# with all pairwise combinations
dfj <- data.frame(t1 = rep(t1, each = length(t2)),
                  t2 = rep(t2, length(t1)))
dfj$z <- dsinvchisq(dfj$t2^2, n-1, s2) * 2*dfj$t2 * dnorm(dfj$t1, my, dfj$t2/sqrt(n))
# breaks for plotting the contours
cl <- seq(1e-5, max(dfj$z), length.out = 6)

















windshieldy_test <- c(13.357, 14.928, 14.896, 14.820)
data = windshieldy_test
n <- length(data)
my <- mean(data)
s2 <- var(data)   # sum((data - my)^2)/(n-1)


tl1 <- c(10, 25)
df1 <- data.frame(t1 = seq(tl1[1], tl1[2], length.out = 1000))



df1$pm_mu <- dt((df1$t1 - my) / sqrt(s2/n), n-1) / sqrt(s2/n)

prob = 0.95
my + qt(c((1-prob)/2, prob+(1-prob)/2), df=n-1, ncp=(1+1/n)*(s2))


df2 <- gather(df1, grp, p, -t1)

p2 <- ggplot(data = df2) +
  geom_line(aes(t1, p, color = grp)) +
  coord_cartesian(xlim = c(-42, 42))
