data {
  int<lower=0> n;             // number of data points
  int<lower=0> k;             // number of covariates
  int<lower=0> y[n];          // time series
  real a1;               // prior mean for the initial state
  real P1;            // prior covariance for the initial state
  real sd_prior_means;   // prior means for the sd parameters
  real sd_prior_sds;     // prior sds for the sd parameters
  vector[k] beta_prior_means; // prior means for the beta parameters
  vector[k] beta_prior_sds;   // prior sds for the sd parameters
  matrix[n, k] xreg;          // covariates
}

parameters {
  real<lower=0> theta;     // sd parameters for level
  vector[k] beta;             // regression coefficients 
  // instead of working directly with true states level and slope
  // it is often suggested use standard normal variables in sampling
  // and reconstruct the true parameters in transformed parameters block
  // this should make sampling more efficient although coding the model 
  // is less intuitive...
  vector[n] level_std;        // N(0, 1) level noise
}

transformed parameters {
  vector[n] xbeta = xreg * beta;
  vector[n] level;
  // construct the actual states
  level[1] = a1 + sqrt(P1) * level_std[1];
  for(t in 2:n) {
    level[t] = level[t-1] + theta * level_std[t];
  }
}

model {
  // priors for theta and beta
  target += normal_lpdf(theta | sd_prior_means, sd_prior_sds);
  target += normal_lpdf(beta | beta_prior_means, beta_prior_sds);
  // standardised noise terms
  target += normal_lpdf(level_std | 0, 1);
  // Poisson likelihood
  target += poisson_log_lpmf(y | xbeta + level);
} 
