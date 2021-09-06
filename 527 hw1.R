
## write simulation code to validate the orthogonal optimization

## Set up true regression model -----
## y = b1*x1 + b2*x2 + epsilon
## epsilon ~ N(0, sigma)
set.seed(2021)
library(mvtnorm)
n = 100
b1 = 2
b2 = 5
sigma = 0.1

n.sim = 10000
b2_list = c()
b2_multi = c()
for ( i in 1:n.sim){
  X = rmvnorm(n, mean=c(0.5,1.5), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
  x1 = X[,1]
  x2 = X[,2]
  epsilon = rnorm(n, 0, sd=sigma)
  
  y = b1*x1 + b2*x2 + epsilon
  ## Obtain b2 via orthogonal optimization-----
  ### 1. regress y ~ x1 -----
  ## exclude components explained by x1 from y
  fit1 <- lm(y ~ x1-1)
  e_yx1 <- fit1$residuals
  
  ### 2. regress x2 ~ x1 -----
  ## exclude components explained by x1 from x2
  fit2 <- lm(x2 ~ x1-1)
  e_x2x1 <- fit2$residuals
  
  ### 3. regress e_yx1 ~ e_x2x1 -----
  fit3 <- lm(e_yx1 ~ e_x2x1-1)
  b2_est <- summary(fit3)$coefficients[1,1]
  b2_list <- c(b2_list,b2_est)
  
}


## estimated b2 is quite close to the true b2:
hist(b2_list)
mean(b2_list) # 5.00008

## averaged biasness of estimated b2 is almost 0:
mean(b2_list-b2) ## 7.954301e-05

## RMSE of estimated b2 is:
sqrt(mean((b2_list-b2)^2)) # 0.007




