data {
  int<lower=0> n;
  int y[n];
  real x1;
  real P1;
  real prior_mean;
  real prior_sd;
}

parameters {
  real<lower=0> sd_x;
  vector[n] x_raw;
}

transformed parameters {
  vector[n] x;
  x[1] = x1 + sqrt(P1) * x_raw[1];
  for(t in 2:n) {
    x[t] = x[t-1] + sd_x * x_raw[t];
  }
}

model {
  target += normal_lpdf(x_raw | 0, 1);
  target += poisson_lpmf(y | exp(x));
}
