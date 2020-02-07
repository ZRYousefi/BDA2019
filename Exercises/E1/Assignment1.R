install.packages("devtools")
devtools::install_github("avehtari/BDA_course_Aalto",
                         subdir = "rpackage")

install.packages("markmyassignment")

library(markmyassignment)
exercise_path <- "https://github.com/avehtari/BDA_course_Aalto/blob/master/exercises/tests/ex1.yml"
set_assignment(exercise_path)

x <- seq(0, 1, length = 101)
mu = 0.2
sigma <- 0.1
a <- mu*((mu*(1-mu))/(sigma^2) - 1)
b <- alpha*(1-mu)/mu
f_1 <- dbeta(x, a, b)
plot(x, fx)


f_2 <- rbeta(1000,a, b) 
hist(f_2)

sample_mean <- mean(f_2)
sample_sigma <- var(f_2)

quantile(f_2)


"==========================="


boxes <- matrix(c(2,4,1,5,1,3), ncol = 2,
                dimnames = list(c("A", "B", "C"), c("red", "white")))
boxes

prob_box = c(0.4, 0.1, 0.5)
p_r_box <- function(boxes){
  p = c()
  for(i in 1:nrow(boxes)){
    p[i] <- boxes[i,1]/sum(boxes[i,])
  }
  return(p) 
}
p_r_box(boxes)

p_red <- function(boxes){
  p = 0
  for(i in 1:nrow(boxes)){
   p <- p + (p_r_box(boxes)[i]*prob_box[i])
  }
  return(p)
}

p_red(boxes)


p_box <- function(boxes){
  p = c()
  for(i in 1:nrow(boxes)){
    p[i] <- (p_r_box(boxes)[i]*prob_box[i])/p_red(boxes)
  }
  return(p)
}

p_box(boxes)


mark_my_assignment()

"==========================="

p_identical_twin <- function(fraternal_prob = 1/125, identical_prob = 1/300){
  p_identical_twin <- 0.5*identical_prob / (0.5*identical_prob + 0.25*fraternal_prob)
  return(p_identical_twin)
}

p_identical_twin(fraternal_prob = 1/125, identical_prob = 1/300)

p_identical_twin(1/150, 1/400)

mark_my_assignment()
