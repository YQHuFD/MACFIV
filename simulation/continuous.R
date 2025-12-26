#' ============================================================================
#' MACFIV Continuous Instrumental Variables Simulation Study
#' ============================================================================
#' 
#' A comprehensive simulation study to evaluate the MACFIV method for 
#' nonlinear causal effect estimation with continuous instrumental variables
#' (genetic instruments with correlated structure).
#'
#' This simulation extends the standard IV simulation to include:
#'   - Continuous instruments with AR(1) correlation structure
#'   - Different functional forms of nonlinear causal effects
#'   - Varying degrees of pleiotropy (direct effects on outcome)
#'   - Performance assessment via bias metrics
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
#'
#' This approach is robust to weak and pleiotropic instruments because:
#' - Only uses information from instrument-exposure relationships
#' - Ignores direct pleiotropic effects on outcome
#' - Combines multiple submodels for stability

MMACF <- function(X, Z) {
  
  n <- nrow(Z)  # Number of observations
  p <- ncol(Z)  # Number of instruments
  
  # ========================================================================
  # STEP 1: Calculate absolute Pearson correlations
  # ========================================================================
  # Measure strength of relationship between exposure and each instrument
  pearson <- rep(0, p)
  
  for (j in 1:p) {
    pearson[j] <- abs(cor(X, Z[, j]))
  }
  
  # ========================================================================
  # STEP 2: Sort instruments by correlation strength (descending)
  # ========================================================================
  # Instruments with stronger correlations are prioritized
  sorted_indices <- order(pearson, decreasing = TRUE)
  
  # ========================================================================
  # STEP 3: Construct nested sequence of submodels
  # ========================================================================
  K_Q <- c(1:p)  # Model sizes: 1, 2, ..., p
  Q <- length(K_Q)  # Total number of submodels
  
  # Store predictions from each submodel
  X_models <- matrix(0, nrow = n, ncol = Q)
  
  # For each model size q, fit model using top q instruments
  for (q in K_Q) {
    Z_q <- as.matrix(Z[, sorted_indices[1:q]])
    # Fit: X ~ Z_q (no intercept)
    X_models[, q] <- predict(lm(X ~ Z_q - 1))
  }
  
  # ========================================================================
  # STEP 4: Estimate error variance (sigma^2)
  # ========================================================================
  # Use full model residuals to estimate error variance
  sigma2 <- sum((residuals(lm(X ~ Z - 1)))^2) / (n - p)
  
  # ========================================================================
  # STEP 5: Define Mallows' criterion objective function
  # ========================================================================
  # Balance model fit against complexity
  objective_function <- function(weights) {
    fitted_values <- X_models %*% weights
    rss <- sum((X - fitted_values)^2)
    # Mallows' criterion: RSS + penalty for model complexity
    mallows <- rss + 2 * sigma2 * weights %*% c(1:Q)
    return(mallows)
  }
  
  # ========================================================================
  # STEP 6: Optimize weights via constrained optimization
  # ========================================================================
  # Constraints:
  #   - sum(weights) = 1  (convex combination)
  #   - weights >= 0      (non-negative)
  
  result <- solnp(
    pars = rep(1 / Q, Q),      # Initial: uniform weights
    fun = objective_function,
    eqfun = function(w) sum(w) - 1,  # Equality: sum = 1
    eqB = 0,
    LB = rep(0, Q),             # Lower bounds: weights >= 0
    control = list(trace = 0)   # Suppress output
  )
  
  # ========================================================================
  # STEP 7: Compute confounding factors as residuals
  # ========================================================================
  weights <- result$pars  # Optimal weights
  X_MMA <- X_models %*% weights  # Weighted model prediction
  CF <- X - X_MMA  # Confounding factors as residuals
  
  return(CF)
}


# ============================================================================
# MAIN FUNCTION: simulate_macfiv_continuous()
# ============================================================================
#' Simulation Study for MACFIV with Continuous Instruments
#'
#' Comprehensive simulation to evaluate MACFIV performance with continuous
#' instrumental variables having AR(1) correlation structure (realistic
#' scenario for genetic instruments with linkage disequilibrium).
#'
#' @param n Integer. Sample size.
#' @param p Integer. Number of instruments.
#' @param s Integer. Number of pleiotropic instruments (direct effect on outcome).
#' @param num_simulations Integer. Number of simulation replicates.
#' @param save_path Character. Directory path for saving results.
#'                   Default: './macfiv_continuous_results/'
#' @param seed Integer. Random seed for reproducibility. Default: NULL.
#'
#' @details
#' Data Generation Process:
#'
#' 1. **Instruments (Z)**:
#'    - Generated from multivariate normal with AR(1) correlation structure
#'    - Correlation(Z_i, Z_j) = 0.5^|i-j| (realistic linkage disequilibrium)
#'    - Standardized to mean 0, variance 1
#'
#' 2. **Unobserved Confounder (U)**:
#'    - U ~ N(0,1) independent of Z
#'
#' 3. **Error Terms**:
#'    - epsilon_X ~ N(0,1) independent error for exposure
#'    - epsilon_Y ~ N(0,1) independent error for outcome
#'    - delta_1 = U + epsilon_X (error for exposure equation)
#'    - delta_2 = delta_1 + epsilon_Y (error for outcome equation)
#'
#' 4. **Exposure (X)**:
#'    - X = Z*gamma + delta_1
#'    - gamma_j = sqrt(2/n) for all j (weak instrument design)
#'
#' 5. **Outcome (Y)**:
#'    - Y = f(X) + Z*alpha + delta_2
#'    - alpha_j = 1 for j = 1,...,s (s pleiotropic instruments)
#'    - alpha_j = 0 for j = s+1,...,p (no pleiotropy)
#'    - f(X) depends on functional form case
#'
#' Functional Forms Tested:
#' 1. null: f(X) = 0, f'(X) = 0
#' 2. linear: f(X) = X, f'(X) = 1
#' 3. quad: f(X) = 0.01*X^2, f'(X) = 0.02*X
#' 4. sin: f(X) = sin(0.1*X), f'(X) = 0.1*cos(0.1*X)
#' 5. exp: f(X) = exp(0.1*X), f'(X) = 0.1*exp(0.1*X)
#' 6. log: f(X) = log(X - min(X) + 1), f'(X) = 1/(X - min(X) + 1)
#' 7. mixed: f(X) = cos(0.1*X) + exp(0.1*X) + 0.01*X^2
#'          f'(X) = -0.1*sin(0.1*X) + 0.1*exp(0.1*X) + 0.02*X
#'
#' MACFIV Methodology:
#' 1. Estimate confounding factors using MMACF function
#' 2. Create B-spline basis to flexibly approximate f(X)
#' 3. Construct design matrix: [B-splines | Instruments | Confounding Factors]
#' 4. Apply SCAD penalty:
#'    - Zero penalty on B-spline terms (preserve smooth function)
#'    - Unit penalty on instruments and confounding factors
#' 5. Select optimal lambda using BIC
#' 6. Extract estimated causal function and marginal effect
#'
#' Output Files:
#' For each functional form, generates CSV file with columns:
#' - bias_original: Bias in function estimate E[f_hat(X) - f(X)]
#' - bias_derivative: Bias in derivative estimate E[f'_hat(X) - f'(X)]
#'
#' @examples
#' \dontrun{
#' # Run simulation with continuous instruments
#' # n=2000, p=100, s=40 (40% pleiotropic), 100 replicates
#' simulate_macfiv_continuous(n = 2000, p = 100, s = 40, 
#'                            num_simulations = 100,
#'                            save_path = './results/', 
#'                            seed = 2024)
#' }
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import ncvreg
#' @import splines
#' @import splines2
#' @export

simulate_macfiv_continuous <- function(n, p, s, num_simulations,
                                        save_path = './macfiv_continuous_results/',
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
  # MAIN SIMULATION LOOP (over functional forms)
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
    # PARALLEL SIMULATION LOOP (over replicates)
    # ====================================================================
    foreach(sim = 1:num_simulations, .combine = 'rbind',
            .packages = c('Rsolnp', 'MASS', 'mgcv', 'ncvreg',
                         'glmnet', 'splines', 'splines2')) %dopar% {
      
      # ==================================================================
      # DATA GENERATION: Continuous Instruments with AR(1) Structure
      # ==================================================================
      
      # Create AR(1) correlation matrix: rho^|i-j|
      # This mimics linkage disequilibrium in genetic data
      mean_vector <- rep(0, p)
      indices <- 1:p
      cov_matrix <- outer(indices, indices, 
                         FUN = function(i, j) 0.5^abs(i - j))
      
      # Generate multivariate normal with AR(1) structure
      Z <- MASS::mvrnorm(n, mu = mean_vector, Sigma = cov_matrix)
      Z <- scale(Z)  # Standardize to mean 0, var 1
      
      # Generate unobserved confounder
      U <- rnorm(n)
      epsilon_X <- rnorm(n)
      epsilon_Y <- rnorm(n)
      delta_1 <- U + epsilon_X
      delta_2 <- delta_1 + epsilon_Y
      
      # Weak instrument design: gamma_j = sqrt(2/n)
      gamma <- rep(sqrt(2/n), p)
      
      # Pleiotropic effects: first s instruments have direct effect on Y
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
      
      # Generate outcome with causal effect and pleiotropy
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
      # This captures the confounding that arises from unobserved confounder
      CF <- MMACF(X, Z)
      
      # Construct design matrix with three components:
      # 1. B-spline basis: captures nonlinear causal effect f(X)
      # 2. Instruments Z: captures direct instrument effects
      # 3. Confounding factors CF: accounts for confounding from pleiotropic effects
      design_matrix <- cbind(B_spline_basis, Z, CF)
      
      # Define penalty structure for SCAD regularization
      lambda_values <- 10^seq(3, -3, -0.1)
      penalty_factors <- rep(0, M + p + 1)
      
      # Zero penalty on B-spline terms (indices 1:M)
      # Preserve smooth function shape by not penalizing
      
      # Unit penalty on instruments (indices M+1:M+p)
      # Allow selection of relevant instruments
      penalty_factors[(M + 1):(M + p)] <- 1
      
      # Zero penalty on confounding factor (index M+p+1)
      # Always include confounding factor
      
      # Fit regularized model with SCAD penalty
      fit <- ncvreg(design_matrix, Y,
                    penalty = "SCAD",
                    penalty.factor = penalty_factors,
                    lambda = lambda_values)
      
      # Select optimal lambda using BIC
      # BIC balances fit quality against model complexity
      bic_values <- BIC(fit)
      best_lambda_index <- which.min(bic_values)
      best_lambda <- fit$lambda[best_lambda_index]
      
      # Fit final model with selected lambda
      model <- ncvfit(design_matrix, Y,
                      penalty = "SCAD",
                      penalty.factor = penalty_factors,
                      lambda = best_lambda)
      
      # ==================================================================
      # EXTRACT AND COMPUTE RESULTS
      # ==================================================================
      
      # Extract B-spline coefficients (first M coefficients)
      spline_coef <- model$beta[1:M]
      
      # Estimate causal function: f(X) = B_spline_basis %*% spline_coef
      f_hat <- B_spline_basis %*% spline_coef
      
      # Estimate marginal effect (derivative): df/dX = B_derivative %*% spline_coef
      derivative_hat <- B_derivative %*% spline_coef
      
      # Calculate bias in function estimate
      bias_original <- mean(f_hat - f_X)
      
      # Calculate bias in derivative estimate
      bias_derivative <- mean(derivative_hat - derivative)
      
      # Store results
      results_macfiv <- matrix(c(bias_original, bias_derivative),
                               nrow = 1)
      
      # Return results to main process
      results_macfiv
      
    }  # End foreach (parallel) loop
    
    # Write results to file
    write.table(results_macfiv, file = filename_macfiv, sep = "\t",
                row.names = FALSE, col.names = FALSE, append = TRUE)
    
    cat("Completed case:", case, "\n")
    
  }  # End case loop
  
  # ========================================================================
  # CLEANUP AND SUMMARY
  # ========================================================================
  stopImplicitCluster()
  stopCluster(cl)
  
  cat("\n", "="*70, "\n")
  cat("Simulation Summary:\n")
  cat("- Sample size (n):", n, "\n")
  cat("- Number of instruments (p):", p, "\n")
  cat("- Number of pleiotropic instruments (s):", s, "\n")
  cat("- Number of replicates:", num_simulations, "\n")
  cat("- Functional forms tested:", length(caseset), "\n")
  cat("- Results saved to:", save_path, "\n")
  cat("="*70, "\n\n")
  
}


# ============================================================================
# EXAMPLE USAGE
# ============================================================================
# Uncomment to run simulation

# Run simulation study with continuous instruments
# Parameters:
#   n=2000 (sample size)
#   p=100 (100 instruments)
#   s=40 (40 pleiotropic instruments, 40% pleiotropy rate)
#   100 simulation replicates
#   Save to './macfiv_continuous_results/'
#   Seed=2024 for reproducibility
#
# simulate_macfiv_continuous(n = 2000, p = 100, s = 40,
#                            num_simulations = 100,
#                            save_path = './macfiv_continuous_results/',
#                            seed = 2024)
