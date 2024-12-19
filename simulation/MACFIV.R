library(Rsolnp)
library(MASS)
library(mgcv)
library(ncvreg)
library(glmnet)
library(splines)
library(splines2)

MMACF <- function(X, Z){
  
  n <- nrow(Z)
  p <- ncol(Z)
  
  pearson <- rep(0,p)
  
  for (j in 1:p){
    pearson[j] <- abs(cor(X,Z[ ,j]))
  }
  
  sorted_indices <- order(pearson, decreasing = TRUE)
  
  K_Q <- c(1:p)
  
  #the number of submodels
  Q <- length(K_Q)     
  
  X_models <- matrix(0, nrow = n, ncol = Q)
  
  for (q in K_Q){
    
    Z_q <- as.matrix(Z[ ,sorted_indices[1:q]])
    X_models[,q] <- predict(lm(X ~ Z_q - 1))
    
  }
  
  sigma2 <- sum((residuals(lm(X ~ Z - 1)))^2)/(n - p)
  
  # Mallows criterion
  objective_function <- function(weights) {
    fitted_values <- X_models %*% weights
    sum((X - fitted_values)^2) + 2 * sigma2 * weights %*% c(1:Q)  # Mallows准则
  }
  
  # restrictive condition
  result <- solnp(pars = rep(1 / Q, Q),  
                  fun = objective_function,
                  eqfun = function(w) sum(w) - 1,  
                  eqB = 0,
                  LB = rep(0, Q),
                  control = list(trace = 0) 
  )
  
  
  #optimal weight
  weights <- result$pars
  X_MMA <- X_models %*% weights
  CF <- X - X_MMA
  
  return(CF)
}








MACFIV <- function(X, Y, Z, degree = 3, M = 5){
  
  n <- nrow(Z)
  p <- ncol(Z)
  
  #B-spline 
  B_spline_basis <- bs(X, degree = degree, df=M, intercept = FALSE)
  
  #B-spline of derivative
  B_derivative <- dbs(X, degree = degree, df=M, intercept = FALSE)
  
  #model-averaged control function
  CF <- MMACF(X,Z)
  
  #hyperparameters of SCAD
  lambda_values <- 10^seq(3, -3, -0.1)
  penalty_factors <- rep(0, M+p+1)
  penalty_factors[(M+1):(M+p)] <- 1
  
  #SCAD
  fit <- ncvreg(cbind(B_spline_basis, Z, CF), Y, penalty="SCAD", penalty.factor=penalty_factors, lambda=lambda_values)
  
  #BIC
  bic_values <- BIC(fit)
  best_lambda_index <- which.min(bic_values)
  best_lambda <- fit$lambda[best_lambda_index]
  
  model <- ncvfit(cbind(B_spline_basis, Z, CF), Y, penalty = "SCAD", penalty.factor = penalty_factors, lambda = best_lambda)
  
  #B-spline fit to g
  f_hat <- B_spline_basis %*% model$beta[1:M]
  
  #B-spline fit to g'
  derivative_hat <- B_derivative %*% model$beta[1:M]
}
