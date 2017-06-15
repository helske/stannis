data {
  int<lower=0> n;             // number of data points
  int<lower=0> y[n];          // time series
  vector[2] a1;               // prior mean for the initial state
  matrix[2, 2] P1;            // prior covariance for the initial state
  vector[2] sd_prior_means;   // prior means for the sd parameters
  vector[2] sd_prior_sds;     // prior sds for the sd parameters
}

parameters {
  real<lower=0> theta[2];     // sd parameters for level and slope
  // instead of working directly with true states level and slope
  // it is often suggested use standard normal variables in sampling
  // and reconstruct the true parameters in transformed parameters block
  // this should make sampling more efficient although coding the model 
  // is less intuitive...
  vector[n] level_std;        // N(0, 1) level noise
  vector[n] slope_std;        // N(0, 1) slope noise
}

transformed parameters {
  vector[n] level;
  vector[n] slope;
  // construct the actual states
  // note that although P1 was allowed to have general form here
  // it is assumed that it is diagonal... laziness (and covers typical cases)
  level[1] = a1[1] + sqrt(P1[1,1]) * level_std[1];
  slope[1] = a1[2] + sqrt(P1[2,2]) * slope_std[1];
  for(t in 2:n) {
    level[t] = level[t-1] + slope[t-1] + theta[1] * level_std[t];
    slope[t] = slope[t-1] + theta[2] * slope_std[t];
  }
}

model {
  // priors for theta
  theta ~ normal(sd_prior_means, sd_prior_sds);
  // standardised noise terms
  level_std ~ normal(0, 1);
  slope_std ~ normal(0, 1);
  // Poisson likelihood
  y ~ poisson_log(level);
} 
