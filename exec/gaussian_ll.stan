data {
  int<lower=0> n;
  real y[n];
  real a1;
  real P1;
  real sd_prior_means[2];
  real sd_prior_sds[2];
}

parameters {
  real<lower=0> theta[2];
  vector[n] x_raw;
}

transformed parameters {
  vector[n] x;
  x[1] = a1 + sqrt(P1) * x_raw[1];
  for(t in 2:n) {
    x[t] = x[t-1] + theta[1] * x_raw[t];
  }
}

model {

  theta ~ normal(sd_prior_means, sd_prior_sds);
  x_raw ~ normal(0, 1);
  y ~ normal(x, theta[2]);
}
