data {
  int<lower=0> n;
  int<lower=0> k;
  int<lower=0> y[n];
  vector[2] x1;
  matrix[2, 2] P1;
  vector[2] sd_prior_means;
  vector[2] sd_prior_sds;
  vector[k] beta_prior_means;
  vector[k] beta_prior_sds;
  matrix[n, k] xreg;
}

parameters {
  real<lower=0> theta[2];
  vector[k] beta;
  vector[n] level;
  vector[n] slope;
}

transformed parameters {
  vector[n] xbeta = xreg * beta;

}
model {
  target += normal_lpdf(theta | sd_prior_means, sd_prior_sds);
  target += normal_lpdf(beta | beta_prior_means, beta_prior_sds);
  target += normal_lpdf(level[1] | x1[1], sqrt(P1[1, 1]));
  target += normal_lpdf(slope[1] | x1[2], sqrt(P1[2, 2]));
  for(t in 2:n) {
    target += normal_lpdf(level[t] | level[t - 1] + slope[t - 1], theta[1]);
    target += normal_lpdf(slope[t] | slope[t - 1], theta[2]);
  }
  target += poisson_log_lpmf(y | xbeta + level);
}
