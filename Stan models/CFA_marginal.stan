// Functions block
functions {
  // Compute sign adjustments for confirmatory models
  vector compute_sign_adjustment(int J, int M, matrix Lambda, matrix L_ind) {
        vector[M] sign_adjustment = rep_vector(1.0, M);
        array[J] int item_used = rep_array(0, J);

        for (m in 1:M) {
            int j_m = 0;
            for (j in 1:J) {
                if (L_ind[j, m] != 0 && item_used[j] == 0) {
                    j_m = j;
                    item_used[j] = 1;
                    break;
                }
            }
            if (j_m > 0 && Lambda[j_m, m] < 0) {
                sign_adjustment[m] = -1.0;
            }
        }
        return sign_adjustment;
    }
}

// Data block
data {
  int<lower=0> N;             // Number of observations
  int<lower=0> J;             // Number of items
  int<lower=0> M;             // Number of latent factors
  matrix[N, J] Y;             // Observed data
  matrix[J, M] L_ind;         // Lambda index matrix
}

// Parameters block
parameters {
  vector[J] mu;                   // Model-implied mean vector
  matrix[J, M] L_unc;             // Lambda matrix
  vector<lower=0>[J] Theta;                // Uniqueness
  cholesky_factor_corr[M] Psi_chol;    // Cholesky decomposition (correlation)
}

// Transformed parameters block
transformed parameters{
  // Constrain Lambda matrix
  matrix[J, M] L_mat = L_unc .* L_ind;
  // Model-implied latent-variable correlation matrix
  corr_matrix[M] Psi_b = multiply_lower_tri_self_transpose(Psi_chol);
  // Model-implied covariance matrix
  cov_matrix[J] Sigma;
  Sigma = L_mat * Psi_b * L_mat' + diag_matrix(Theta);
}

// Model block
model {
  // Model priors
  mu ~ normal(0, 10);
  to_vector(L_unc) ~ normal(0, 5);
  Theta ~ gamma(1, 0.5);
  Psi_chol ~ lkj_corr_cholesky(1);
  
  // Model log-likelihood
  for(i in 1:N){
    Y[i,] ~ multi_normal(mu, Sigma);
  }
}

// Generated quantities block
generated quantities {
  // Compute sign adjustment
  vector[M] sign_adjustment = compute_sign_adjustment(J, M, L_mat, L_ind);
  
  // Adjusted factor loadings
  matrix[J, M] Lambda = L_mat;
  for (m in 1:M) { 
    Lambda[, m] *= sign_adjustment[m]; 
    }
  
  // Adjusted latent factor correlation matrix
  corr_matrix[M] Psi = Psi_b;
  for(m1 in 1:M){
    for(m2 in 1:M){ 
      Psi[m1, m2] *= sign_adjustment[m1] * sign_adjustment[m2]; 
      }
    }

  // Model log-likelihood
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = multi_normal_lpdf(Y[i,] | mu, Sigma);
  }
}
