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
  theta ~ normal(sd_prior_means, sd_prior_sds);
  x_raw ~ normal(0, 1);
  y ~ poisson_log(x);
}
