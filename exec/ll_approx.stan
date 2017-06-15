functions {
  #include "common_functions.stan"
}

data {
  int<lower=0> n;
  vector[n] y;
  vector[1] a1;
  matrix[1, 1] P1;
  real sd_prior_means;
  real sd_prior_sds;
  vector[n+1] initial_mode;
  int distribution;
  int max_iter;
  real conv_tol;
}

transformed data {

  int<lower=0> m = 1;
  row_vector[m] Zt;
  matrix[m, m] Tt;
  vector[n] xbeta = rep_vector(0.0, n);
  Zt[1] = 1.0;
  Tt[1, 1] = 1.0;
}

parameters {
  real<lower=0> theta;
}

transformed parameters {
  vector[3 * n + 1] approx_results;
  matrix[m, m] Rt;
  Rt[1, 1] = theta^2;
  approx_results = approx(y, a1, P1, Zt, Tt, Rt, initial_mode, xbeta, distribution,
  max_iter, conv_tol);
}

model {
  theta ~ normal(sd_prior_means, sd_prior_sds);
  target += approx_results[3 * n + 1];
}

generated quantities {
  real jacobian;
  jacobian = -log(theta);
}

