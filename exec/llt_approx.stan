functions {
  #include "common_functions.stan"
}

data {
  int<lower=0> n;
  vector[n] y;
  vector[2] x1;
  matrix[2, 2] P1;
  vector[2] sd_prior_means;
  vector[2] sd_prior_sds;
  vector[n+1] initial_mode;
  int distribution;
}

transformed data {

  int<lower=0> m = 2;
  row_vector[m] Zt;
  matrix[m, m] Tt;
  vector[n] xbeta = rep_vector(0.0, n);
  Zt[1] = 1.0;
  Zt[2] = 0.0;
  Tt[1, 1] = 1.0;
  Tt[2, 1] = 0.0;
  Tt[1, 2] = 1.0;
  Tt[2, 2] = 1.0;
}


parameters {
  real<lower=0> theta[2];
}

transformed parameters {
  vector[3 * n + 1] approx_results;
  matrix[m, m] Rt = rep_matrix(0.0, m, m);
  Rt[1, 1] = theta[1]^2;
  Rt[2, 2] = theta[2]^2;
  approx_results = approx(y, x1, P1, Zt, Tt, Rt, initial_mode, xbeta, distribution);
}
model {
  target += normal_lpdf(theta | sd_prior_means, sd_prior_sds);
  target += approx_results[3 * n + 1];
}

generated quantities {
  real jacobian;
  jacobian = -sum(log(theta));
}
