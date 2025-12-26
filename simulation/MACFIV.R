#' ============================================================================
#' MACFIV: A Novel Framework for Nonlinear Causal Inference with
#'         Many Weak and Pleiotropic Genetic Instruments
#' ============================================================================
#'
#' This package implements the MACFIV method for estimating nonlinear causal
#' effects in the presence of many weak and pleiotropic genetic instruments.
#' It combines multiple adaptive confounding factor estimation with 
#' semiparametric regression using B-splines and sparse variable selection.
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
# FUNCTION 1: normalize() - Data Normalization
# ============================================================================
#' Normalize Data Matrix
#'
#' Normalize columns of a data matrix using min-max scaling based on a base matrix.
#'
#' @param X Numeric matrix to be normalized
#' @param base Numeric matrix used as reference for normalization statistics.
#'             If NULL, X is used as reference (default).
#'
#' @return Normalized matrix with same dimensions as X

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
      # Min-max normalization: (x - mean(base)) / (max(base) - min(base))
      X[, ip] <- (Xcol - mean(base)) / (max(base) - min(base))
    }
  }
  
  return(X)
}


# ============================================================================
# FUNCTION 2: MMACF() - Multiple Adaptive Confounding Factors Estimation
# ============================================================================
#' Estimate Multiple Adaptive Confounding Factors
#'
#' This function estimates confounding factors by constructing multiple
#' submodels based on instrument strength (correlation) and optimally
#' combining them using Mallows' criterion.
#'
#' @param X Numeric vector. The exposure variable (n x 1).
#' @param Z Numeric matrix. The instrumental variables (n x p).
#'
#' @details
#' The MMACF algorithm:
#' 1. Calculates Pearson correlations between exposure X and each instrument
#' 2. Orders instruments by decreasing correlation strength
#' 3. Constructs Q nested submodels where model q uses top q instruments
#' 4. Uses Mallows' criterion to optimally weight these submodels
#' 5. Returns residuals as estimated confounding factors
#'
#' The Mallows' criterion balances model fit (RSS) against complexity:
#'   C_p = RSS + 2*sigma2*p
#' where p is the number of parameters.
#'
#' @return
#' Numeric vector of confounding factors (residuals from optimal weighted model).
#' These factors capture the confounding variation not explained by instruments.
#'
#' @references
#' Mallows, C. L. (1973). "Some comments on Cp". Technometrics, 15(4), 661-675.

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
  # Instruments with stronger correlations are used first
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
  objective_function <- function(weights) {
    # Weighted combination of all submodel predictions
    fitted_values <- X_models %*% weights
    
    # Residual sum of squares
    rss <- sum((X - fitted_values)^2)
    
    # Mallows' criterion: RSS + penalty for model complexity
    # sum(weights) gives the effective number of "models" selected
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
    LB = rep(0, Q),             # Lower bounds
    control = list(trace = 0)   # Suppress output
  )
  
  # ========================================================================
  # STEP 7: Compute final confounding factors
  # ========================================================================
  weights <- result$pars  # Optimal weights
  X_MMA <- X_models %*% weights  # Weighted model prediction
  CF <- X - X_MMA  # Confounding factors as residuals
  
  return(CF)
}


# ============================================================================
# FUNCTION 3: MACFIV() - Complete Nonlinear Causal Effect Estimation
# ============================================================================
#' MACFIV: Nonlinear Causal Effect Estimation with Multiple Adaptive
#'         Confounding Factors and Semiparametric Regression
#'
#' Estimates the nonlinear causal effect of exposure X on outcome Y,
#' accounting for many weak and pleiotropic genetic instruments.
#'
#' @param X Numeric vector. The exposure variable (n x 1).
#' @param Y Numeric vector. The outcome variable (n x 1).
#' @param Z Numeric matrix. The instrumental variables (n x p).
#' @param degree Integer. Degree of B-spline basis (default: 3, cubic).
#' @param df Integer. Degrees of freedom for B-spline basis (default: 5).
#'
#' @return
#' A list containing:
#'   \item{f_hat}{Estimated causal function values}
#'   \item{derivative_hat}{Estimated causal derivative (marginal effect)}
#'   \item{spline_coef}{B-spline coefficients}
#'   \item{lambda}{Selected tuning parameter (lambda)}
#'   \item{model}{Fitted ncvreg model object}
#'
#' @details
#' The MACFIV framework consists of 5 main steps:
#'
#' 1. **Estimate Confounding Factors (MMACF)**
#'    - Construct multiple instrument submodels by correlation strength
#'    - Use Mallows' criterion to optimally weight submodels
#'    - Extract confounding factors from residuals
#'
#' 2. **Construct B-spline Basis**
#'    - Create smooth basis functions to flexibly approximate nonlinear relationship
#'    - Use cubic B-splines with 5 degrees of freedom
#'
#' 3. **Combine Design Matrix**
#'    - Include: B-spline basis + all instruments + confounding factors
#'    - This accounts for direct effects, instrument effects, and confounding
#'
#' 4. **Apply Sparse Variable Selection (SCAD)**
#'    - Use SCAD penalty to select relevant features
#'    - Set zero penalty on B-spline coefficients (preserve smooth function)
#'    - Set unit penalty on instruments and confounding factors
#'    - Choose lambda by BIC
#'
#' 5. **Extract Causal Function and Derivative**
#'    - Predicted values from B-spline terms = estimated causal function
#'    - Derivative of B-spline basis = marginal effect of exposure
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' n <- 500
#' p <- 50
#' Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- Z[, 1:10] %*% rnorm(10) + rnorm(n)
#' Y <- sin(X/5) + Z[, 1:5] %*% rnorm(5) + rnorm(n)
#'
#' # Estimate causal effect
#' result <- MACFIV(X, Y, Z)
#' 
#' # Extract results
#' f_hat <- result$f_hat
#' derivative_hat <- result$derivative_hat
#' }
#'
#' @import mgcv
#' @import ncvreg
#' @import splines
#' @import splines2
#' @export

MACFIV <- function(X, Y, Z, degree = 3, df = 5) {
  
  n <- nrow(Z)  # Number of observations
  p <- ncol(Z)  # Number of instruments
  
  # ========================================================================
  # STEP 1: ESTIMATE CONFOUNDING FACTORS using MMACF
  # ========================================================================
  # Construct multiple adaptive confounding factors from instruments
  CF <- MMACF(X, Z)
  
  # ========================================================================
  # STEP 2: CONSTRUCT B-SPLINE BASIS FUNCTIONS
  # ========================================================================
  # Create smooth basis functions for flexible nonlinear approximation
  # bs() creates the basis functions
  # dbs() creates their derivatives (for marginal effects)
  
  B_spline_basis <- bs(X, degree = degree, df = df, intercept = FALSE)
  B_derivative <- dbs(X, degree = degree, df = df, intercept = FALSE)
  
  M <- ncol(B_spline_basis)  # Number of B-spline basis functions
  
  # ========================================================================
  # STEP 3: CONSTRUCT DESIGN MATRIX
  # ========================================================================
  # Combine three components:
  #   1. B-spline basis: captures nonlinear causal effect f(X)
  #   2. Instruments Z: captures direct instrument effects
  #   3. Confounding factors CF: accounts for confounding from pleiotropic effects
  
  design_matrix <- cbind(B_spline_basis, Z, CF)
  
  # ========================================================================
  # STEP 4: SPARSE VARIABLE SELECTION with SCAD PENALTY
  # ========================================================================
  # Define penalty structure:
  #   - 0 penalty on B-spline terms (preserve smooth function shape)
  #   - 1 penalty on instruments and confounding factors (select relevant ones)
  
  lambda_values <- 10^seq(3, -3, -0.1)  # Range of lambda values to search
  
  penalty_factors <- rep(0, M + p + 1)
  # Zero penalty for B-spline terms (indices 1:M)
  # Unit penalty for instruments and confounding factor (indices M+1:M+p+1)
  penalty_factors[(M + 1):(M + p + 1)] <- 1
  
  # Fit regularized model with SCAD penalty across lambda sequence
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
  
  # ========================================================================
  # STEP 5: EXTRACT CAUSAL FUNCTION AND MARGINAL EFFECT
  # ========================================================================
  # Extract B-spline coefficients (indices 1:M)
  spline_coef <- model$beta[1:M]
  
  # Estimated causal function: f(X) = B_spline_basis * spline_coef
  f_hat <- B_spline_basis %*% spline_coef
  
  # Estimated marginal effect (derivative): df/dX = B_derivative * spline_coef
  derivative_hat <- B_derivative %*% spline_coef
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  return(list(
    f_hat = as.numeric(f_hat),
    derivative_hat = as.numeric(derivative_hat),
    spline_coef = spline_coef,
    lambda = best_lambda,
    model = model,
    confounding_factors = CF
  ))
}


# ============================================================================
# EXAMPLE USAGE
# ============================================================================
# Uncomment to run example simulation
#
# # Generate synthetic data
# set.seed(123)
# n <- 500
# p <- 100
# 
# # Generate instruments
# Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
# Z <- scale(Z)  # Standardize
# 
# # Generate exposure with instrument effects
# gamma <- rep(sqrt(2/n), p)
# X <- Z %*% gamma + rnorm(n)
# 
# # Generate outcome with nonlinear causal effect
# f_X <- sin(0.1 * X)
# Y <- f_X + Z[, 1:5] %*% rnorm(5) + rnorm(n)
# 
# # Apply MACFIV
# result <- MACFIV(X, Y, Z, degree = 3, df = 5)
# 
# # Extract results
# f_hat <- result$f_hat
# derivative_hat <- result$derivative_hat
# 
# # Compute estimation error
# mae_f <- mean(abs(f_X - f_hat))
# mae_deriv <- mean(abs(0.1*cos(0.1*X) - derivative_hat))
# 
# cat("MAE for function estimate:", mae_f, "\n")
# cat("MAE for derivative estimate:", mae_deriv, "\n")
