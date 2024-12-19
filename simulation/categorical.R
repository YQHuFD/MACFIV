library(parallel)
library(foreach)
library(doParallel)

library(MASS)
library(mgcv)
library(ncvreg)
library(glmnet)
library(splines)
library(splines2)


simulation_invalid <- function(n, p, s, num_simulations){
  
  save_path <- '/vhome/shidapeng/cd/MACFIR/categorical/(1000,100,10)/'
  num_cores <- 40
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("normalize", "MMACF"))
  
  methods <- c('TSP', 'TSP-SCAD', 'PolyMR', 'CF', 'MACFIR')
  caseset <- c('null', 'linear', 'quad', 'sin', 'exp', 'log', 'mixed')
  
  for (case in caseset){
    
    filename_TSP <- paste0(save_path, case, " (", n, ",", p, ",", s, ") ", "TSP", ".csv")
    filename_TSP_SCAD <- paste0(save_path, case, " (", n, ",", p, ",", s, ") ", "TSP-SCAD", ".csv")
    filename_PolyMR <- paste0(save_path, case, " (", n, ",", p, ",", s, ") ", "PolyMR", ".csv")
    filename_CF <- paste0(save_path, case, " (", n, ",", p, ",", s, ") ", "CF", ".csv")
    filename_MACFIR <- paste0(save_path, case, " (", n, ",", p, ",", s, ") ", "MACFIR", ".csv")
    
    if (!file.create(filename_TSP)) {
      stop(paste("Error: Unable to create file", filename_TSP))
    }
    if (!file.create(filename_TSP_SCAD)) {
      stop(paste("Error: Unable to create file", filename_TSP_SCAD))
    }
    if (!file.create(filename_PolyMR)) {
      stop(paste("Error: Unable to create file", filename_PolyMR))
    }
    if (!file.create(filename_CF)) {
      stop(paste("Error: Unable to create file", filename_CF))
    }
    if (!file.create(filename_MACFIR)) {
      stop(paste("Error: Unable to create file", filename_MACFIR))
    }
    
    
    titles <- matrix(c("mae_original", "rmse_original", "mae_derivative", "rmse_derivative"), nrow = 1)
    write.table(titles, file = filename_TSP, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(titles, file = filename_TSP_SCAD, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(titles, file = filename_PolyMR, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(titles, file = filename_CF, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(titles, file = filename_MACFIR, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    
    
    
    foreach (sim = 1:num_simulations, .combine = 'rbind', .packages = c('Rsolnp', 'MASS', 'mgcv', 'ncvreg', 'glmnet', 'splines', 'splines2', 'PolyMR')) %dopar% {
      
      tau <- matrix(rbinom(n * p, size = 1, prob = 0.3), nrow = n, ncol = p) 
      xi <- matrix(rbinom(n * p, size = 1, prob = 0.3), nrow = n, ncol = p) 
      Z <- tau + xi

      U <- rnorm(n)  # Unobserved confounder
      epsilon_X <- rnorm(n)  # Error term for X
      epsilon_Y <- rnorm(n)  # Error term for Y
      delta_1 <- U + epsilon_X
      delta_2 <- delta_1 + epsilon_Y
      
      gamma <- rep(sqrt(2/n),p)

      alpha <- rep(0, p)
      if (s != 0){
        alpha[1:s] <- 1
      }
      X <- Z %*% gamma + delta_1
      
      
      if (case == 'null'){
        f_X <- 0
        derivative <- rep(0,n)
      }
      
      if (case == 'linear'){
        f_X <- X
        derivative <- rep(1,n)
      }
      
      if (case == 'quad'){
        f_X <- 0.1*X^2
        derivative <- 0.2*X
      }
      
      if (case == 'sin'){
        f_X <- sin(X)
        derivative <- cos(X)
      }
      
      if (case == 'exp'){
        f_X <- exp(0.5*X)
        derivative <- 0.5*exp(0.5*X)
      }
      
      if (case == 'log'){
        f_X <- log((X-min(X))+1)
        derivative <- 1/((X-min(X))+1)
      }
      
      if (case == 'mixed'){
        f_X <- cos(X) + exp(0.5*X) + 0.1*X^2
        derivative <- -sin(X) + 0.5*exp(0.5*X) +  0.2*X
      }
      
      Y <- f_X + Z %*% alpha + delta_2
      
      n <- nrow(Z)
      p <- ncol(Z)

      degree <- 3 
      M  <- 5   
      B_spline_basis <- bs(X, degree = degree, df=M, intercept = FALSE)
      B_derivative <- dbs(X, degree = degree, df=M, intercept = FALSE)
      
      
      for (method in methods){
        
        if (method == 'TSP'){
          
          hat_B <- matrix(0, n, M)
          for (i in 1:M) {
            fit <- lm(B_spline_basis[, i] ~ Z)
            hat_B[, i] <- predict(fit)
          }
          D <- hat_B
          coef <- solve(t(D) %*% D) %*% t(D) %*% Y
          
          f_hat <- D %*% coef
          mae_original_TSP <- mean(abs(f_X - f_hat))
          rmse_original_TSP <- sqrt(mean((f_X - f_hat)^2))
          
          derivative_hat <- B_derivative %*% coef
          mae_derivative_TSP <- mean(abs(derivative - derivative_hat))
          rmse_derivative_TSP <- sqrt(mean((derivative - derivative_hat)^2))
          
          results_TSP <- matrix(c(mae_original_TSP, rmse_original_TSP, mae_derivative_TSP, rmse_derivative_TSP), nrow = 1)
        }
        
        
        if (method == 'TSP-SCAD'){
          
          lambda_values <- 10^seq(3, -3, -0.1)
          penalty_factors <- rep(0, M+p)
          penalty_factors[(M+1):(M+p)] <- 1
          
          hat_B <- matrix(0, n, M)
          for (i in 1:M) {
            fit <- lm(B_spline_basis[, i] ~ Z)
            hat_B[, i] <- predict(fit)
          }
          fit <- ncvreg(cbind(hat_B, Z), Y, penalty="SCAD", penalty.factor=penalty_factors, lambda=lambda_values)
          bic_values <- BIC(fit)
          best_lambda_index <- which.min(bic_values)
          best_lambda <- fit$lambda[best_lambda_index]
          
          model <- ncvfit(cbind(hat_B, Z), Y, penalty = "SCAD", penalty.factor = penalty_factors, lambda = best_lambda)
          
          f_hat <- hat_B %*% model$beta[1:M]
          mae_original_TSP_SCAD <- mean(abs(f_X - f_hat))
          rmse_original_TSP_SCAD <- sqrt(mean((f_X - f_hat)^2))
          
          derivative_hat <- B_derivative %*% model$beta[1:M]
          mae_derivative_TSP_SCAD <- mean(abs(derivative - derivative_hat))
          rmse_derivative_TSP_SCAD <- sqrt(mean((derivative - derivative_hat)^2))
          
          results_TSP_SCAD <- matrix(c(mae_original_TSP_SCAD, rmse_original_TSP_SCAD, mae_derivative_TSP_SCAD, rmse_derivative_TSP_SCAD), nrow = 1)
        }
        
        
        if (method == 'PolyMR'){
          
          fit <- polymr(X,Y,Z)
          coefficients <- unname(fit$polymr$outcome_model$coefficients)[1:11]
          phenotypes_summary <- fit$phenotypes_summary
          mean_X <- phenotypes_summary$mean[phenotypes_summary$phenotype == "Exposure"]
          sd_X <- phenotypes_summary$sd[phenotypes_summary$phenotype == "Exposure"]
          mean_Y <- phenotypes_summary$mean[phenotypes_summary$phenotype == "Outcome"]
          sd_Y <- phenotypes_summary$sd[phenotypes_summary$phenotype == "Outcome"]
  
          f_hat <- 0
          derivative_hat <- 0
          
          for (i in 1:length(coefficients)) {
            f_hat <- f_hat + coefficients[i] * X_scale^(i - 1)
          }
          
          
          for (i in 2:length(coefficients)) {
            derivative_hat <- derivative_hat + coefficients[i] * (i - 1) * X_scale^(i - 2)
          }
          
          mae_original_polymr <- mean(abs(f_X - f_hat))
          rmse_original_polymr <- sqrt(mean((f_X - f_hat)^2))
          
          mae_derivative_polymr <- mean(abs(derivative - derivative_hat))
          rmse_derivative_polymr <- sqrt(mean((derivative - derivative_hat)^2))
          
          results_polymr <- matrix(c(mae_original_polymr, rmse_original_polymr, mae_derivative_polymr, rmse_derivative_polymr), nrow = 1)
        }
        
        
        if (method == 'CF'){
          
          CF <- as.numeric(lm(X ~ Z)$residuals)
          D <- cbind(B_spline_basis, CF)
          coef <- solve(t(D) %*% D) %*% t(D) %*% Y
          
          f_hat <- B_spline_basis %*% coef[1:M]
          mae_original_CF <- mean(abs(f_X - f_hat))
          rmse_original_CF <- sqrt(mean((f_X - f_hat)^2))
          
          derivative_hat <- B_derivative %*% coef[1:M]
          mae_derivative_CF <- mean(abs(derivative - derivative_hat))
          rmse_derivative_CF <- sqrt(mean((derivative - derivative_hat)^2))
          
          results_CF <- matrix(c(mae_original_CF, rmse_original_CF, mae_derivative_CF, rmse_derivative_CF), nrow = 1)
        }
        
        
        if (method == 'MACFIR'){
          
          lambda_values <- 10^seq(3, -3, -0.1)
          CF <- MMACF(X,Z)
          
          p <- ncol(Z)
          penalty_factors <- rep(0, M+p+1)
          penalty_factors[(M+1):(M+p)] <- 1
          fit <- ncvreg(cbind(B_spline_basis, Z, CF), Y, penalty="SCAD", penalty.factor=penalty_factors, lambda=lambda_values)
          bic_values <- BIC(fit)
          best_lambda_index <- which.min(bic_values)
          best_lambda <- fit$lambda[best_lambda_index]
          
          model <- ncvfit(cbind(B_spline_basis, Z, CF), Y, penalty = "SCAD", penalty.factor = penalty_factors, lambda = best_lambda)
          
          f_hat <- B_spline_basis %*% model$beta[1:M]
          mae_original_MACFIR <- mean(abs(f_X - f_hat))
          rmse_original_MACFIR <- sqrt(mean((f_X - f_hat)^2))
          
          derivative_hat <- B_derivative %*% model$beta[1:M]
          mae_derivative_MACFIR <- mean(abs(derivative - derivative_hat))
          rmse_derivative_MACFIR <- sqrt(mean((derivative - derivative_hat)^2))
          
          results_MACFIR <- matrix(c(mae_original_MACFIR, rmse_original_MACFIR, mae_derivative_MACFIR, rmse_derivative_MACFIR), nrow = 1)
        }
      }
      
      write.table(results_TSP, file = filename_TSP, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(results_TSP_SCAD, file = filename_TSP_SCAD, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(results_polymr, file = filename_PolyMR, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(results_CF, file = filename_CF, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(results_MACFIR, file = filename_MACFIR, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
  
  stopImplicitCluster()
  stopCluster(cl)
}

