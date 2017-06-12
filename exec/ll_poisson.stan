data {
  int<lower=0> n;
  int y[n];
  real a1;
  real P1;
  real sd_prior_means;
  real sd_prior_sds;
}

parameters {
  real<lower=0> theta;
  vector[n] x_raw;
}

transformed parameters {
  vector[n] x;
  x[1] = a1 + sqrt(P1) * x_raw[1];
  for(t in 2:n) {
    x[t] = x[t-1] + theta * x_raw[t];
  }
}

model {
  target += normal_lpdf(theta | sd_prior_means, sd_prior_sds);
  target += normal_lpdf(x_raw | 0, 1);
  target += poisson_log_lpmf(y | x);
}
