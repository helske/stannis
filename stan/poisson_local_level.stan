data {
  int<lower=0> n;
  int y[n];
  real x1;
  real P1;
}

parameters {
  real<lower=0> sd_x;
  vector[n] x;
}

model {
  target += normal_lpdf(sd_x | 0, 10);
  target += normal_lpdf(x[1] | x1, P1);
  for(t in 2:n) {
    target += normal_lpdf(x[t] | x[t-1], sd_x);
  }
  target += poisson_lpmf(y | exp(x));
}
