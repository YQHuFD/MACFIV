#' ============================================================================
#' MACFIV Simulation Study
#' ============================================================================
#' 
#' A comprehensive simulation study to evaluate the MACFIV method for 
#' nonlinear causal effect estimation with many weak and pleiotropic 
#' genetic instruments.
#'
#' The simulation examines the performance under various scenarios:
#'   - Different functional forms: null, linear, quadratic, sine, exponential, 
#'                                 logarithmic, and mixed functions
#'   - Different sample sizes (n), number of instruments (p), and 
#'     number of pleiotropic instruments (s)
#'   - Measures: bias in function estimation and derivative (marginal effect)
#'
#' ============================================================================

# Required packages
library(Rsolnp)
library(MASS)
library(mgcv)
library(ncvreg)
library(glmnet)
library(splines)
library(splines2)
library(parallel)
library(foreach)
library(doParallel)


# ============================================================================
# HELPER FUNCTION: normalize()
# ============================================================================
#' Normalize Data Matrix
#'
#' Normalize columns of a data matrix using min-max scaling.
#'
#' @param X Numeric matrix to be normalized.
#' @param base Numeric matrix used as reference for normalization statistics.
#'             If NULL (default), X is used as reference.
#'
#' @return Normalized matrix with same dimensions as X.
#'
#' @details
#' Min-max normalization formula:
#'   X_normalized[,i] = (X[,i] - mean(base[,i])) / (max(base[,i]) - min(base[,i]))

normalize <- function(X, base = NULL) {
  X <- as.matrix(X)
  
  if (is.null(base)) {
    base <- X
  }
  base <- as.matrix(base)
  
  p <- dim(X)[2]
  
  # Normalize each column using min-max scaling
  for (ip in 1:p) {
    Xcol <- X[, ip]
    base.col <- base[, ip]
    
    if (var(Xcol) > 0) {
      X[, ip] <- (Xcol - mean(base)) / (max(base) - min(base))
    }
  }
  
  return(X)
}


# ============================================================================
# FUNCTION: MMACF() - Multiple Adaptive Confounding Factors
# ============================================================================
#' Estimate Multiple Adaptive Confounding Factors
#'
#' Constructs multiple instrument submodels and optimally combines them
#' using Mallows' criterion to estimate confounding factors.
#'
#' @param X Numeric vector. The exposure variable (n x 1).
#' @param Z Numeric matrix. The instrumental variables (n x p).
#'
#' @return Numeric vector of confounding factors (residuals).
#'
#' @details
#' Algorithm steps:
#' 1. Calculate Pearson correlations between exposure and each instrument
#' 2. Order instruments by decreasing correlation strength
#' 3. Construct Q nested submodels where model q uses top q instruments
#' 4. Use Mallows' criterion to optimally weight these submodels
#' 5. Return residuals as estimated confounding factors
#'
#' Mallows' criterion: C_p = RSS + 2*sigma2*p

MMACF <- function(X, Z) {
  
  n <- nrow(Z)  # Number of observations
  p <- ncol(Z)  # Number of instruments
  
  # ========================================================================
  # STEP 1: Calculate absolute Pearson correlations
  # ========================================================================
  pearson <- rep(0, p)
  
  for (j in 1:p) {
    pearson[j] <- abs(cor(X, Z[, j]))
  }
  
  # ========================================================================
  # STEP 2: Sort instruments by correlation strength (descending)
  # ========================================================================
  sorted_indices <- order(pearson, decreasing = TRUE)
  
  # ========================================================================
  # STEP 3: Construct nested sequence of submodels
  # ========================================================================
  K_Q <- c(1:p)  # Model sizes: 1, 2, ..., p
  Q <- length(K_Q)  # Total number of submodels
  
  # Store predictions from each submodel
  X_models <- matrix(0, nrow = n, ncol = Q)
  
  for (q in K_Q) {
    Z_q <- as.matrix(Z[, sorted_indices[1:q]])
    # Fit: X ~ Z_q (no intercept)
    X_models[, q] <- predict(lm(X ~ Z_q - 1))
  }
  
  # ========================================================================
  # STEP 4: Estimate error variance (sigma^2)
  # ========================================================================
  sigma2 <- sum((residuals(lm(X ~ Z - 1)))^2) / (n - p)
  
  # ========================================================================
  # STEP 5: Define Mallows' criterion objective function
  # ========================================================================
  objective_function <- function(weights) {
    fitted_values <- X_models %*% weights
    rss <- sum((X - fitted_values)^2)
    mallows <- rss + 2 * sigma2 * weights %*% c(1:Q)
    return(mallows)
  }
  
  # ========================================================================
  # STEP 6: Optimize weights via constrained optimization
  # ========================================================================
  result <- solnp(
    pars = rep(1 / Q, Q),
    fun = objective_function,
    eqfun = function(w) sum(w) - 1,  # Constraint: sum(weights) = 1
    eqB = 0,
    LB = rep(0, Q),  # Constraint: weights >= 0
    control = list(trace = 0)
  )
  
  # ========================================================================
  # STEP 7: Compute confounding factors
  # ========================================================================
  weights <- result$pars
  X_MMA <- X_models %*% weights
  CF <- X - X_MMA
  
  return(CF)
}


# ============================================================================
# MAIN FUNCTION: simulate_macfiv()
# ============================================================================
#' Simulation Study for MACFIV Method
#'
#' Comprehensive simulation to evaluate MACFIV performance under various
#' scenarios with different functional forms and confounding structures.
#'
#' @param n Integer. Sample size.
#' @param p Integer. Number of instruments.
#' @param s Integer. Number of pleiotropic instruments (direct effect on outcome).
#' @param num_simulations Integer. Number of simulation replicates.
#' @param save_path Character. Directory path for saving results.
#'                   Default: './macfiv_results/'
#' @param seed Integer. Random seed for reproducibility. Default: NULL.
#'
#' @details
#' Data Generation Process:
#' - Instruments Z: Binary matrix from Bernoulli(0.3) + Bernoulli(0.3)
#' - Exposure X: X = Z*gamma + delta_1, where gamma = sqrt(2/n)
#' - Unobserved confounder U and errors independent of Z
#' - Outcome Y: Y = f(X) + Z*alpha + delta_2
#'   where alpha[1:s] = 1 (s pleiotropic instruments)
#'
#' Functional Forms Tested:
#' 1. null: f(X) = 0
#' 2. linear: f(X) = X
#' 3. quad: f(X) = 0.01*X^2
#' 4. sin: f(X) = sin(0.1*X)
#' 5. exp: f(X) = exp(0.1*X)
#' 6. log: f(X) = log(X - min(X) + 1)
#' 7. mixed: f(X) = cos(0.1*X) + exp(0.1*X) + 0.01*X^2
#'
#' Output Files:
#' For each functional form, generates CSV file with columns:
#' - bias_original: Bias in function estimate E[f_hat(X) - f(X)]
#' - bias_derivative: Bias in derivative estimate E[f'_hat(X) - f'(X)]
#'
#' @examples
#' \dontrun{
#' # Run simulation with n=500, p=100, s=10, 100 replicates
#' simulate_macfiv(n = 500, p = 100, s = 10, num_simulations = 100,
#'                 save_path = './results/', seed = 2024)
#' }
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import ncvreg
#' @import splines
#' @import splines2
#' @export

simulate_macfiv <- function(n, p, s, num_simulations, 
                            save_path = './macfiv_results/', 
                            seed = NULL) {
  
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create directory if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # ========================================================================
  # PARALLEL COMPUTING SETUP
  # ========================================================================
  num_cores <- detectCores() - 1  # Use all available cores except one
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export functions to cluster nodes
  clusterExport(cl, c("normalize", "MMACF"), envir = environment())
  
  # ========================================================================
  # SIMULATION PARAMETERS
  # ========================================================================
  # Different functional forms to test
  caseset <- c('null', 'linear', 'quad', 'sin', 'exp', 'log', 'mixed')
  
  # B-spline settings
  degree <- 3   # Cubic splines
  M <- 5        # Number of B-spline basis functions
  
  # ========================================================================
  # MAIN SIMULATION LOOP
  # ========================================================================
  for (case in caseset) {
    
    # Construct output filename
    filename_macfiv <- paste0(
      save_path, case, " (", n, ",", p, ",", s, ") MACFIV.csv"
    )
    
    # Create output file
    if (!file.create(filename_macfiv)) {
      stop(paste("Error: Unable to create file", filename_macfiv))
    }
    
    # Write header
    titles <- matrix(c("bias_original", "bias_derivative"), nrow = 1)
    write.table(titles, file = filename_macfiv, sep = "\t",
                row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # ====================================================================
    # PARALLEL SIMULATION LOOP
    # ====================================================================
    foreach(sim = 1:num_simulations, .combine = 'rbind',
            .packages = c('Rsolnp', 'MASS', 'mgcv', 'ncvreg', 
                         'glmnet', 'splines', 'splines2')) %dopar% {
      
      # ==================================================================
      # DATA GENERATION
      # ==================================================================
      
      # Generate binary instruments from Bernoulli mixture
      tau <- matrix(rbinom(n * p, size = 1, prob = 0.3), 
                    nrow = n, ncol = p)
      xi <- matrix(rbinom(n * p, size = 1, prob = 0.3), 
                   nrow = n, ncol = p)
      Z <- tau + xi
      Z <- scale(Z)  # Standardize
      
      # Generate unobserved confounder and error terms
      U <- rnorm(n)
      epsilon_X <- rnorm(n)
      epsilon_Y <- rnorm(n)
      delta_1 <- U + epsilon_X
      delta_2 <- delta_1 + epsilon_Y
      
      # Instrument effects on exposure
      gamma <- rep(sqrt(2/n), p)
      
      # Pleiotropic effects (direct effects on outcome)
      alpha <- rep(0, p)
      if (s != 0) {
        alpha[1:s] <- 1
      }
      
      # Generate exposure
      X <- Z %*% gamma + delta_1
      
      # ==================================================================
      # GENERATE OUTCOME BASED ON FUNCTIONAL FORM
      # ==================================================================
      
      if (case == 'null') {
        f_X <- 0
        derivative <- rep(0, n)
      }
      
      if (case == 'linear') {
        f_X <- X
        derivative <- rep(1, n)
      }
      
      if (case == 'quad') {
        f_X <- 0.01 * X^2
        derivative <- 0.02 * X
      }
      
      if (case == 'sin') {
        f_X <- sin(0.1 * X)
        derivative <- 0.1 * cos(0.1 * X)
      }
      
      if (case == 'exp') {
        f_X <- exp(0.1 * X)
        derivative <- 0.1 * exp(0.1 * X)
      }
      
      if (case == 'log') {
        f_X <- log((X - min(X)) + 1)
        derivative <- 1 / ((X - min(X)) + 1)
      }
      
      if (case == 'mixed') {
        f_X <- cos(0.1 * X) + exp(0.1 * X) + 0.01 * X^2
        derivative <- -0.1 * sin(0.1 * X) + 0.1 * exp(0.1 * X) + 0.02 * X
      }
      
      # Generate outcome
      Y <- f_X + Z %*% alpha + delta_2
      
      # ==================================================================
      # MACFIV METHOD IMPLEMENTATION
      # ==================================================================
      
      # Create B-spline basis functions and derivatives
      B_spline_basis <- bs(X, degree = degree, df = M, 
                           intercept = FALSE)
      B_derivative <- dbs(X, degree = degree, df = M, 
                          intercept = FALSE)
      
      # Estimate confounding factors using MMACF
      CF <- MMACF(X, Z)
      
      # Construct design matrix: [B-splines | Instruments | Confounding Factors]
      design_matrix <- cbind(B_spline_basis, Z, CF)
      
      # Define penalty structure
      lambda_values <- 10^seq(3, -3, -0.1)
      penalty_factors <- rep(0, M + p + 1)
      # Zero penalty on B-spline terms (preserve smooth function)
      # Unit penalty on instruments and confounding factor (select relevant ones)
      penalty_factors[(M + 1):(M + p)] <- 1
      
      # Fit regularized model with SCAD penalty
      fit <- ncvreg(design_matrix, Y, 
                    penalty = "SCAD",
                    penalty.factor = penalty_factors,
                    lambda = lambda_values)
      
      # Select optimal lambda using BIC
      bic_values <- BIC(fit)
      best_lambda_index <- which.min(bic_values)
      best_lambda <- fit$lambda[best_lambda_index]
      
      # Fit final model with selected lambda
      model <- ncvfit(design_matrix, Y,
                      penalty = "SCAD",
                      penalty.factor = penalty_factors,
                      lambda = best_lambda)
      
      # ==================================================================
      # EXTRACT RESULTS
      # ==================================================================
      
      # Extract B-spline coefficients
      spline_coef <- model$beta[1:M]
      
      # Estimate causal function
      f_hat <- B_spline_basis %*% spline_coef
      
      # Estimate marginal effect (derivative)
      derivative_hat <- B_derivative %*% spline_coef
      
      # Calculate bias
      bias_original <- mean(f_hat - f_X)
      bias_derivative <- mean(derivative_hat - derivative)
      
      # Store results
      results_macfiv <- matrix(c(bias_original, bias_derivative), 
                               nrow = 1)
      
      # Return results to main process
      results_macfiv
      
    }  # End foreach loop
    
    # Write results to file
    write.table(results_macfiv, file = filename_macfiv, sep = "\t",
                row.names = FALSE, col.names = FALSE, append = TRUE)
    
    cat("Completed case:", case, "\n")
    
  }  # End case loop
  
  # ========================================================================
  # CLEANUP
  # ========================================================================
  stopImplicitCluster()
  stopCluster(cl)
  
  cat("\nSimulation completed. Results saved to:", save_path, "\n")
  
}


# ============================================================================
# EXAMPLE USAGE
# ============================================================================
# Uncomment to run simulation

# Run simulation study
# Parameters: n=2000, p=100, s=10 (10% pleiotropic instruments), 
#             100 replicates, save to './macfiv_results/'
# simulate_macfiv(n = 2000, p = 100, s = 10, num_simulations = 100,
#                 save_path = './macfiv_results/', seed = 2024)
