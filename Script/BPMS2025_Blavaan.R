# ==============================================================================
# HEADER: Script Information
# ==============================================================================
# Script Name: BPM_Blavaan.R
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
library(blavaan)       # Estimate bayesian models
library(ggplot2)       # Density plots
library(bayesplot)     # Plot results
library(cowplot)       # Arrange plots

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
Y <- as.data.frame(mvrnorm(n = 300, mu = Nu, Sigma = Sigma))

# Check simulated values
data.frame(sample = colMeans(Y), population = Nu)
data.frame(sample = apply(Y, 2, sd), population = item_sds)
data.frame(sample = cor(Y)[lower.tri(cor(Y))], population = Rho[lower.tri(Rho)])

# Parameters in raw scale
Lambda_unstd <- Lambda * item_sds
uniqueness_unstd <- uniqueness * item_sds^2
# ==============================================================================

# ==============================================================================
# SECTION 4: Estimate bayesian psychometric model with Blavaan
# Parameter estimates, convergence, fit measures and model comparison
# ==============================================================================

# Model syntax: Null model
model.null <- '
V1 ~~ V1 
V2 ~~ V2 
V3 ~~ V3
V4 ~~ V4
V5 ~~ V5
V6 ~~ V6 '

# Model syntax: one latent common-factor
model.1f <- '
F1 =~ V1 + V2 + V3 + V4 + V5 + V6
'
# Model syntax: Two correlated common factors
model.2f <- '
F1 =~ V1 + V2 + V3
F2 =~ V4 + V5 + V6
'

# ----------------------------------------------- #
#    Prior Predictive Checking: default priors    #
# ----------------------------------------------- #

# By-default prior distributions
dpriors()

# Fit an "empty" model with prisamp = TRUE! (ignore all the warnings)
prior_predictive <- bcfa(model.2f, data = Y, std.lv = TRUE, 
                         meanstructure = TRUE, test = "none", 
                         sample = 2e4, prisamp = TRUE, 
                         bcontrol = list(cores = 3)) 

# Save prior distribution
prior.draws <- blavInspect(prior_predictive, "mcmc")

# Prior predictive distribution: item intercepts
color_scheme_set("orange")
mcmc_hist(x = prior.draws, regex_pars = "~1", alpha = .5)

# Prior predictive distribution: factor loadings
color_scheme_set("blue")
mcmc_hist(x = prior.draws, regex_pars = "=~", alpha = .5)

# Prior predictive distribution: latent correlation matrix
color_scheme_set("teal")
mcmc_hist(x = prior.draws, pars = "F1~~F2", alpha = .5)

# Prior predictive distribution: residual variances
color_scheme_set("red")
mcmc_hist(x = prior.draws, pars = paste0("V", 1:6, "~~", "V", 1:6), alpha = .5)

# -------------------------------------------------- #
#    Prior Predictive Checking: specifying priors    #
# -------------------------------------------------- #

# Now, let's modify the prior parameters. For example, Lambda and nu.
(new_priors <- dpriors(nu = "normal(0,15)", lambda = "normal(0,5)", 
                       rho = "beta(2,2)", psi = "gamma(1,.2)"))

# Fit an "empty" model with prisamp = TRUE! (ignore all the warnings)
new_prior_predictive <- bcfa(model.2f, data = Y, std.lv = TRUE, 
                             meanstructure = TRUE, test = "none", 
                             sample = 2e4, prisamp = TRUE, 
                             dp = new_priors, # Added new priors
                             bcontrol = list(cores = 3)) 

# Check the new prior predictive distribution
# Save prior distribution
new.prior.draws <- blavInspect(new_prior_predictive, "mcmc")

# Prior predictive distribution: item intercepts
color_scheme_set("orange")
mcmc_hist(x = new.prior.draws, regex_pars = "~1", alpha = .5)

# Prior predictive distribution: factor loadings
color_scheme_set("blue")
mcmc_hist(x = new.prior.draws, regex_pars = "=~", alpha = .5)

# Prior predictive distribution: latent correlation matrix
color_scheme_set("teal")
mcmc_hist(x = new.prior.draws, pars = "F1~~F2", alpha = .5)

# Prior predictive distribution: residual variances
color_scheme_set("red")
mcmc_hist(x = new.prior.draws, pars = paste0("V", 1:6, "~~", "V", 1:6), alpha = .5)

# ------------------------------------------------ #
#    Prior Predictive Checking: simulating data    #
# ------------------------------------------------ #

# Just simulate 100 datasets 
prior_simdata <- sampleData(new_prior_predictive, nrep = 100, simplify = TRUE)

# Observed vs prior predicted distribution for one item
prior_plots_list <- vector(mode = "list", length = ncol(Y))
for(i in 1:ncol(Y)){
  prior_plots_list[[i]] <- ppc_dens_overlay(
    y = Y[,i], yrep = t(sapply(prior_simdata, function(x) x[,i]))
  )
}

# All posterior plots
plot_grid(plotlist = prior_plots_list)

# ------------------------------- #
#    Bayesian Model Estimation    #
# ------------------------------- #

# Estimate all the models
blavaan.null.fit <- bcfa(model.null, data = Y, burnin = 500, sample = 1000, 
                         meanstructure = TRUE,std.lv = TRUE, bcontrol = list(cores = 3))
blavaan.1f.fit <- bcfa(model.1f, data = Y, burnin = 500, sample = 1000, 
                       meanstructure = TRUE, std.lv = TRUE, bcontrol = list(cores = 3))
blavaan.2f.fit <- bcfa(model.2f, data = Y, burnin = 500, sample = 1000,  
                       meanstructure = TRUE, std.lv = TRUE, bcontrol = list(cores = 3))

# --------------------------------------------------- #
#    Bayesian Convergence and Efficiency Assesment    #
# --------------------------------------------------- #

# Convergence: Potential Scale Reduction Factor (PSRF) must be below 1.05
which(blavInspect(blavaan.1f.fit, what = "rhat") > 1.05)
which(blavInspect(blavaan.2f.fit, what = "rhat") > 1.05)

# Efficiency: Effective Sample Sizes (at least 100 ESS per chain)
min(blavInspect(blavaan.1f.fit, "neff"))
min(blavInspect(blavaan.2f.fit, "neff"))

# Traceplots for parameters
posterior.1f <- blavInspect(blavaan.1f.fit, "mcmc")
posterior.2f <- blavInspect(blavaan.2f.fit, "mcmc")

color_scheme_set("mix-blue-red")
mcmc_trace(x = posterior.1f, regex_pars = "F1=~")
mcmc_trace(x = posterior.2f, regex_pars = c("F1=~", "F2=~"))

# If one chain goes wrong, you can highlight it
color_scheme_set("viridis")
mcmc_trace_highlight(x = posterior.2f, regex_pars = "F1~~F2", highlight = 3)

# --------------------------------------- #
#    Bayesian Approximate Fit Measures    #
# --------------------------------------- #

# Bayesian fit indices: One factor
(blav_fit_1f <- blavFitIndices(object = blavaan.1f.fit, baseline.model = blavaan.null.fit))

# Posterior distribution of bayesian fit indices
color_scheme_set("brightblue")
mcmc_hist(data.frame(blav_fit_1f@indices), alpha = 0.5,
          pars = c("BRMSEA", "BCFI", "BTLI", "BGammaHat"))

# Bayesian fit indices: Two factor
(blav_fit_2f <- blavFitIndices(object = blavaan.2f.fit, baseline.model = blavaan.null.fit))

# Posterior distribution of bayesian fit indices
color_scheme_set("green")
mcmc_hist(data.frame(blav_fit_2f@indices), alpha = 0.5,
          pars = c("BRMSEA", "BCFI", "BTLI", "BGammaHat"))

# ------------------------------- #
#    Bayesian Model Comparison    #
# ------------------------------- #

# Bayesian model comparison with several statistics
blav_com_1f_2f <- blavCompare(object1 = blavaan.1f.fit, object2 = blavaan.2f.fit)

# Watanabe-Akaike Information Criteria
blav_com_1f_2f$waic[[1]]    # Object 1 (i.e., blavaan.1f.fit) WAIC
blav_com_1f_2f$waic[[2]]    # Object 2 (i.e., blavaan.2f.fit) WAIC

# WAIC comparison: best model in the first row (here, model 2)
blav_com_1f_2f$diff_waic

# Leave-One-Out Cross-validation (better than WAIC)
blav_com_1f_2f$loo[[1]]    # Object 1 (i.e., blavaan.1f.fit) WAIC
blav_com_1f_2f$loo[[2]]    # Object 2 (i.e., blavaan.2f.fit) WAIC

# WAIC comparison: best model in the first row (here, model 2)
blav_com_1f_2f$diff_loo

# Log-bayes factor via Laplace-Metropolis Approximation
# Positive values favor object 1 (i.e., blavaan.1f.fit)
blav_com_1f_2f$bf # mll are marginal log-likelihoods

# Just exponentiate it see BF12
exp(blav_com_1f_2f$bf)[1] # Close-to zero.

# Now, BF21 is the inverse of BF12
1/exp(blav_com_1f_2f$bf)[1] 

# ------------------------------ #
#    Plots and model summaries   #
# ------------------------------ #

# Model 1 summary
summary(blavaan.1f.fit, standardized = TRUE, rsquare = TRUE)

# Model 2 summary
summary(blavaan.2f.fit, standardized = TRUE, rsquare = TRUE)

# Posterior distribution: histograms
color_scheme_set("orange")
mcmc_hist(x = posterior.2f, regex_pars = "~1", alpha = .5)

# Posterior distribution: density plots
color_scheme_set("teal")
mcmc_dens(x = posterior.2f, pars = "F1~~F2", alpha = .5)

# Posterior distribution: pairs plots
color_scheme_set("gray")
mcmc_pairs(posterior.2f, regex_pars = "F1=~", diag_fun = "hist", 
           off_diag_fun = "scatter", off_diag_args = list(alpha = 0.5))

# Posterior distribuition: uncertainty intervals
color_scheme_set("pink")
mcmc_intervals(posterior.2f, pars = paste0("V", 1:6, "~~", "V", 1:6))

# Posterior distribuition: uncertainty intervals with areas
color_scheme_set("blue")
mcmc_areas(posterior.2f, regex_pars = c("F1=~", "F2=~"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "median"
)

# --------------------------------------- #
#    Posterior Predictive Model Ckecks    #
# --------------------------------------- #

# I'ts a bit slow, but it works!
# Estimated time: 2 or 3 minutes per posterior_check
posterior_checks_1f <- ppmc(object = blavaan.1f.fit, thin = 1, fit.measures = c("srmr","chisq"))
posterior_checks_2f <- ppmc(object = blavaan.2f.fit, thin = 1, fit.measures = c("srmr","chisq"))

# PPMC: Likelihood-Ratio test (i.e. chisq)
plot(posterior_checks_2f, element = "chisq")
plot(posterior_checks_1f, element = "chisq")
hist(posterior_checks_2f, element = "chisq")
hist(posterior_checks_1f, element = "chisq")

# PPMC: Standardized Root Mean Residuals
plot(posterior_checks_2f, element = "srmr")
plot(posterior_checks_1f, element = "srmr")
hist(posterior_checks_2f, element = "srmr")
hist(posterior_checks_1f, element = "srmr")

# Changing "discrepancy function": model-implied correlation-residual matrix
res_cors <- function(fit){ lavaan::lavResiduals(fit, zstat = FALSE, summary = FALSE)$cov } 
posterior_cor_residuals_1f <- ppmc(object = blavaan.1f.fit, thin = 1, discFUN = res_cors)
posterior_cor_residuals_2f <- ppmc(object = blavaan.2f.fit, thin = 1, discFUN = res_cors)

# PPMC: All correlations
unique_cors <- cbind(paste0("V", combn(1:6,2)[1,]), paste0("V", combn(1:6,2)[2,]))

# One-factor model
par(mfrow=c(3,5))
for(i in 1:nrow(unique_cors)){
  plot(posterior_cor_residuals_1f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
for(i in 1:nrow(unique_cors)){
  hist(posterior_cor_residuals_1f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
par(mfrow=c(1,1))

# Two-factor model
par(mfrow=c(3,5))
for(i in 1:nrow(unique_cors)){
  plot(posterior_cor_residuals_2f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
for(i in 1:nrow(unique_cors)){
  hist(posterior_cor_residuals_2f, element = c(unique_cors[i,1], unique_cors[i,2]))
}
par(mfrow=c(1,1))

# ---------------------------------------------------- #
#    Posterior Predictive Model Ckecks: Reliability    #
# ---------------------------------------------------- #

# Reliability based on composite scores 
bcfa.reliability <- function(fit){ semTools::compRelSEM(fit) }

# Compote PPMC omega statistic
posterior_omega_1f <- ppmc(object = blavaan.1f.fit, thin = 1, discFUN = bcfa.reliability)
posterior_omega_2f <- ppmc(object = blavaan.2f.fit, thin = 1, discFUN = bcfa.reliability)

# Plot reliability PPMC: one common factor
plot(posterior_omega_1f, element = "F1")
hist(posterior_omega_1f, element = "F1")

# Plot reliability PPMC: Two latent factors
par(mfrow=c(1,2))
plot(posterior_omega_2f, element = "F1")
plot(posterior_omega_2f, element = "F2")
par(mfrow=c(1,1))

# Plot reliability PPMC: Two latent factors
par(mfrow=c(1,2))
hist(posterior_omega_2f, element = "F1")
hist(posterior_omega_2f, element = "F2")
par(mfrow=c(1,1))

# Posterior reliability distribution: one common factor
hist(unlist(posterior_omega_1f@obsDist$discFUN1), breaks = 17, 
     xlab = "Reliability", ylab = "Density", 
     main = "Posterior distribution: reliability", 
     col = "powderblue")

# Mean and standard deviation
mean(unlist(posterior_omega_1f@obsDist$discFUN1))
sd(unlist(posterior_omega_1f@obsDist$discFUN1))

# Posterior reliability distribution: Two common factors
reliab.2f <- as.data.frame(do.call(rbind, posterior_omega_2f@obsDist$discFUN1))

# Posterior reliability distribution: one common factor
par(mfrow=c(1,2))
hist(reliab.2f$F1, breaks = 17, 
     xlab = "Reliability: F1", ylab = "Density", 
     main = "Posterior distribution: reliability", 
     col = "peachpuff")
hist(reliab.2f$F2, breaks = 17, 
     xlab = "Reliability: F2", ylab = " ", 
     main = " ", 
     col = "lightblue1")
par(mfrow=c(1,1))

# Mean and standard deviation
colMeans(reliab.2f)
apply(reliab.2f, 2, sd)

# -------------------------------------------------- #
#    Posterior Predictive Ckecks: simulating data    #
# -------------------------------------------------- #

# Just simulate 100 datasets 
posterior_simdata <- sampleData(blavaan.2f.fit, nrep = 100, simplify = TRUE)

# Observed vs prior predicted distribution for one item
posterior_plots_list <- vector(mode = "list", length = ncol(Y))
for(i in 1:ncol(Y)){
  posterior_plots_list[[i]] <- ppc_dens_overlay(
    y = Y[,i], yrep = t(sapply(posterior_simdata, function(x) x[,i]))
    )
}

# All posterior plots
plot_grid(plotlist = posterior_plots_list)

# ==============================================================================
