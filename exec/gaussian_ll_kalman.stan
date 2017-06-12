functions {
 // univariate Kalman filter, return log-likelihood
real gaussian_filter(vector y, vector a1, matrix P1, real Ht,
  row_vector Zt, matrix Tt, matrix Rt) {

  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;

  vector[m] x = a1;
  matrix[m, m] P = P1;

  for (t in 1:n) {
    real F = quad_form(P, Zt') + Ht;
    if (F > 1.0e-8) {
      real v = y[t] - dot_product(Zt, x);
      vector[m] K = P * Zt' / F;
      x = Tt * (x + K * v);
      P = quad_form_sym(P - K * K' * F, Tt') + Rt;
      loglik = loglik - 0.5 * (log(F) + v * v / F);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + Rt;
    }
  }
   return loglik;
  }

}

data {

  int<lower=0> n;
  vector[n] y;
  vector[1] a1;
  matrix[1,1] P1;
  real sd_prior_means[2];
  real sd_prior_sds[2];

}
transformed data {
  int<lower=0> m = 1;
  row_vector[m] Zt;
  matrix[m, m] Tt;
  Zt[1] = 1.0;
  Tt[1, 1] = 1.0;
}
parameters {
  real<lower=0> theta[2];
}

transformed parameters {
  matrix[m, m] Rt;
  Rt[1, 1] = theta[1]^2;
}

model {
  target += normal_lpdf(theta | sd_prior_means, sd_prior_sds);
  target += gaussian_filter(y, a1, P1, theta[2]^2, Zt, Tt, Rt);
}
