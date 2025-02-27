# ==============================================================================
# HEADER: Script Information
# ==============================================================================
# Script Name: BPM_Stan.R
# Author: Ricardo Rey-Sáez
# Position: PhD Student in Psychology, Autonomous University of Madrid, Spain.
# Date Created: 2025-02-27
# Email: ricardoreysaez95@gmail.com
# ==============================================================================

# ==============================================================================
# SECTION 1: Load packages
# Libraries necessary for the script to function
# ==============================================================================

library(MASS)          # Multivariate gaussian distribution
library(cmdstanr)      # Estimate bayesian models in Stan (better than rstan)
library(posterior)     # Extract summaries
library(loo)           # WAIC and LOO estimates

# ==============================================================================

# ==============================================================================
# SECTION 2: Simulate data from CFA model with marginal likelihood approach
# True model: two common latent factors
# ==============================================================================

# Set seed for reproducible simulation
set.seed(2025)

# Population item intercepts
Nu <- rep(c(-1, 0, 1), times = 2)

# Population factor loadings
Lambda <- matrix(c(rep(c(.7, 0), each = 3), 
                   rep(c(0, .7), each = 3)), ncol = 2)

# Population latent factors correlation matrix
Phi <- matrix(c(1, .5, .5, 1), ncol = 2)

# Population item standard deviation (non-standardized)
item_sds <- rep(c(2, 2.5, 3), times = 2)

# Model-implied correlation matrix
communality <- Lambda %*% Phi %*% t(Lambda)
uniqueness <- 1 - diag(communality)
Rho <- communality + diag(uniqueness)

# Population covariance matrix
Sigma <- diag(item_sds) %*% Rho %*% diag(item_sds)

# Simulate data
Y <- mvrnorm(n = 300, mu = Nu, Sigma = Sigma)

# Check simulated values
data.frame(sample = colMeans(Y), population = Nu)
data.frame(sample = apply(Y, 2, sd), population = item_sds)
data.frame(sample = cor(Y)[lower.tri(cor(Y))], population = Rho[lower.tri(Rho)])

# ==============================================================================

# ==============================================================================
# SECTION 3: Estimate bayesian psychometric model with JAGS
# One latent common-factor model
# ==============================================================================

# Load our Stan model
BCFA <- cmdstan_model(stan_file = "Stan models/CFA_marginal.stan")
# Prepare Stan input data
# (this is what this model needs)
# data {
#   int<lower=0> N;             // Number of observations
#   int<lower=0> P;             // Number of items
#   int<lower=0> M;             // Number of latent factors
#   array[N] row_vector [P] Y;  // Observed data
#   matrix[P, M] L_ind;         // Lambda index matrix
# }

# Stan data: 1 factor model (wrong model)
sdata.1f <- list(
  N = 300,
  J = 6,
  M = 1,
  Y = Y,
  L_ind = matrix(rep(1, 6),ncol=1)
)

# Fit the model
BCFA_1f <- BCFA$sample(data = sdata.1f,          # Stan data
                       chains = 4,            # Number of chains
                       parallel_chains = 4,   # Number of parallel chains
                       iter_warmup = 500,     # Adaptation iterations
                       iter_sampling = 1500,  # Sampled iterations
                       adapt_delta = .8,      # 
                       refresh = 500,         # Progress bar
                       init = 0)              # All starting values = 0

# ------------------------ #
#   Parameter estimates    #
# ------------------------ #

# Let's check parameter estimates!
BCFA_1f$summary("Lambda")                # Estimated factor loadings
BCFA_1f$summary("mu")                        # Model-implied mean vector
print(BCFA_1f$summary("Sigma"), n = 1e3)     # Model-implied covariance matrix

# ==============================================================================

# ==============================================================================
# SECTION 4: Estimate bayesian psychometric model with Stan
# Two latent common-factor model
# ==============================================================================

# Stan data: 2 factor model (true model)
sdata.2f <- list(
  N = 300,
  J = 6,
  M = 2,
  Y = Y,
  L_ind = ifelse(Lambda == 0, 0, 1)
)

# Fit the model
BCFA_2f <- BCFA$sample(data = sdata.2f,       # Stan data
                       chains = 4,            # Number of chains
                       parallel_chains = 4,   # Number of parallel chains
                       iter_warmup = 500,     # Adaptation iterations
                       iter_sampling = 1500,  # Sampled iterations
                       refresh = 500,         # Progress bar
                       seed = 2025,           # Reproducible results
                       init = 0)              # All starting values = 0

# ------------------------ #
#   Parameter estimates    #
# ------------------------ #

# Let's check parameter estimates!
BCFA_2f$summary("Lambda")                    # Estimated factor loadings
BCFA_2f$summary("mu")                        # Model-implied mean vector
print(BCFA_2f$summary("Sigma"), n = 1e3)     # Model-implied covariance matrix

# ==============================================================================

# ==============================================================================
# SECTION 5: Model comparison via incremental fit indices
# WAIC and LOO variants
# ==============================================================================

# ------------------------------------- #
#   Bayesian Incremental Fit Indices    #
# ------------------------------------- #

# Watanabee-Akaike Information Criteria
(waic_1f <- loo::waic(x = BCFA_1f$draws("log_lik")))
(waic_2f <- loo::waic(x = BCFA_2f$draws("log_lik")))

# Leave-One-Out Cross-validation
(loo_1f  <- loo::loo(x = BCFA_1f$draws("log_lik")))
(loo_2f  <- loo::loo(x = BCFA_2f$draws("log_lik")))

# Direct loo comparison
loo::loo_compare(loo_1f, loo_2f)

# Pseudo-Bayesian model averaging 
loo_model_weights(list(mod.1f = loo_1f, mod.2f = loo_2f), 
                  method = "pseudobma")

# Pseudo-Bayesian model averaging with bayesian bootstrap
loo_model_weights(list(mod.1f = loo_1f, mod.2f = loo_2f), 
                  method = "pseudobma", BB = TRUE)

# Bayesian model stacking
loo_model_weights(list(mod.1f = loo_1f, mod.2f = loo_2f), 
                  method = "stacking")

# ==============================================================================

