functions {
  #include "common_functions.stan"
}

data {
  int<lower=0> n;
  int<lower=0> k;
  int<lower=0> period;
  vector[n] y;
  vector[2 + period] a1;
  matrix[2 + period, 2 + period] P1;
  vector[3] sd_prior_means;
  vector[3] sd_prior_sds;
  vector[k] beta_prior_means;
  vector[k] beta_prior_sds;
  vector[n+1] initial_mode;
  matrix[n, k] xreg;
  int distribution;
  int max_iter;
  real conv_tol;
}

transformed data {

  int<lower=0> m = 2 + period;
  row_vector[m] Zt = rep_row_vector(0.0, m);
  matrix[m, m] Tt = rep_matrix(0.0, m, m);
  Zt[1] = 1.0;
  Zt[3] = 1.0;
  Tt[1, 1] = 1.0;
  Tt[1, 2] = 1.0;
  Tt[2, 2] = 1.0;
  Tt[3, 3:m] = rep_row_vector(-1.0, m - 2);
  for (i in 4:m) {
    Tt[i, i - 1] = 1.0;
  }
}


parameters {
  real<lower=0> theta[3];
  vector[k] beta;
}

transformed parameters {
  vector[n] xbeta = xreg * beta;
  vector[3 * n + 1] approx_results;
  matrix[m, m] Rt = rep_matrix(0.0, m, m);
  Rt[1, 1] = theta[1]^2;
  Rt[2, 2] = theta[2]^2;
  Rt[3, 3] = theta[3]^2;
  approx_results = approx(y, a1, P1, Zt, Tt, Rt, initial_mode, xbeta, distribution,
  max_iter, conv_tol);
}
model {
  beta ~ normal(beta_prior_means, beta_prior_sds);
  theta ~ normal(sd_prior_means, sd_prior_sds);
  target += approx_results[3 * n + 1];
}

generated quantities {
  real jacobian;
  jacobian = -sum(log(theta));
}
