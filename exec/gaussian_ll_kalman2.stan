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
matrix gaussian_smoother(vector y, vector a1, matrix P1, real var_y,
  row_vector Zt, matrix Tt, matrix Rt) {

  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;

  for (t in 1:n) {
    F[t] = quad_form(P, Zt') + var_y;
    if (F[t] > 1.0e-8) {
      v[t] = y[t] - dot_product(Zt, x);
      K[, t] = P * Zt' / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt') + Rt;
    //Tt * (P - K[,t] * K[,t]' * F[t]) * Tt' + Rt;
      loglik = loglik - 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + Rt;
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[,t+1];
    if(F[t] > 1.0e-8) {
      r[,t] =  Zt' * v[t] / F[t] + (Tt - Tt * K[,t] * Zt)' * tmp;
    } else {
      r[,t] = Tt' * tmp;
    }
  }

  tmpr = r[,1];
  r[,1] = a1 + P1 * tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = Tt * tmp + Rt * tmp2;
  }
  return r[1:m, 1:n];
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
  matrix[m, m] P_L;
  P_L = cholesky_decompose(P1);
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

generated quantities{
    vector[n] y_sim;
    matrix[m, n] a_sim;
    matrix[m, m] R_L;
    vector[m] tmpm;
    R_L = cholesky_decompose(Rt);
    a_sim[1:m, 1] = multi_normal_cholesky_rng(a1, P_L);
    for (t in 1:(n - 1)) {
      tmpm = a_sim[1:m, t];
      a_sim[1:m, t+1]  = multi_normal_cholesky_rng(Tt * tmpm, R_L);
    } 
    for(t in 1:n) {
      y_sim[t] = y[t] - Zt * a_sim[1:m, t] + normal_rng(0, theta[2]);
    }
    a_sim = a_sim + gaussian_smoother(y_sim, a1, P1, theta[2]^2, Zt, Tt, Rt);
}
