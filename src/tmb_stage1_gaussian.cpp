//' ---
//' title: "TMB Stage 1 Gaussian Model"
//' description: "C++ template for Template Model Builder (TMB) to compute the negative log-likelihood of a linear mixed model with random intercept and random slope."
//' ---

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()() {
  using namespace density;

  // --- Data Inputs ---
  DATA_VECTOR(y);             // Level-1 outcome variable (concatenated)
  DATA_VECTOR(z);             // Level-1 predictor variable (concatenated)
  DATA_IVECTOR(cluster_start); // 0-indexed start position of each cluster
  DATA_IVECTOR(cluster_size);  // Number of observations in each cluster

  // --- Parameters ---
  PARAMETER(beta0);  // Fixed intercept
  PARAMETER(beta1);  // Fixed slope
  PARAMETER(var0);   // Random intercept variance
  PARAMETER(cov01);  // Random intercept/slope covariance
  PARAMETER(var1);   // Random slope variance
  PARAMETER(sigma2); // Residual variance

  // Domain checks: Variances must be strictly positive
  if (var0 <= Type(0) || var1 <= Type(0) || sigma2 <= Type(0)) {
    return Type(1e20); // Return large penalty
  }

  // Domain checks: The random-effect covariance matrix G must be positive definite
  Type det_g = var0 * var1 - cov01 * cov01;
  if (det_g <= Type(0)) {
    return Type(1e20); // Return large penalty
  }

  int n_cluster = cluster_start.size();
  
  // Pack fixed effects into a vector
  vector<Type> beta(2);
  beta(0) = beta0;
  beta(1) = beta1;

  // Construct the random-effect covariance matrix G
  matrix<Type> G(2, 2);
  G(0, 0) = var0;
  G(0, 1) = cov01;
  G(1, 0) = cov01;
  G(1, 1) = var1;

  vector<Type> cluster_nll(n_cluster);
  Type nll = Type(0);

  // Loop over each cluster to compute its contribution to the log-likelihood
  for (int i = 0; i < n_cluster; ++i) {
    int start = cluster_start(i);
    int n_i = cluster_size(i);

    // Construct design matrix X and outcome vector y for cluster i
    matrix<Type> X(n_i, 2);
    vector<Type> y_i(n_i);

    for (int j = 0; j < n_i; ++j) {
      X(j, 0) = Type(1); // Intercept column
      X(j, 1) = z(start + j);
      y_i(j) = y(start + j);
    }

    // Marginal covariance for cluster i: V_i = X_i * G * X_i' + sigma^2 * I
    matrix<Type> V_i = X * G * X.transpose();
    for (int j = 0; j < n_i; ++j) {
      V_i(j, j) += sigma2;
    }

    // Compute residuals
    vector<Type> resid_i = y_i - X * beta;
    
    // Evaluate multivariate normal density
    Type contrib = MVNORM(V_i)(resid_i);
    
    cluster_nll(i) = contrib;
    nll += contrib;
  }

  // Report cluster-level contributions (useful for sandwich variance estimator)
  REPORT(cluster_nll);
  
  return nll;
}

