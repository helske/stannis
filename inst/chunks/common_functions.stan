// univariate Kalman filter, return log-likelihood
real gaussian_filter(vector y, vector a1, matrix P1, vector var_y,
  row_vector Zt, matrix Tt, matrix Rt) {

  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;

  vector[m] x = a1;
  matrix[m, m] P = P1;

  for (t in 1:n) {
    real F = quad_form(P, Zt') + var_y[t];
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

// univariate Kalman smoothing
// returns a vector of smoothed means  Z * alpha_t
// and the log-likelihood as a last element of the vector
vector gaussian_smoother(vector y, vector a1, matrix P1, vector var_y,
  row_vector Zt, matrix Tt, matrix Rt) {

  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;
  vector[n+1] mode;
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;

  for (t in 1:n) {
    F[t] = quad_form(P, Zt') + var_y[t];
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

  for(t in 1:n) {
    mode[t] = Zt * r[,t];
  }
  mode[n+1] = loglik;

  return mode;
}

// Gaussian approximation of Poisson model
vector approx(vector y, vector a1, matrix P1, row_vector Zt,
  matrix Tt, matrix Rt, vector mode_, vector xbeta, int distribution, 
  int max_iter, real conv_tol) {

  int n = rows(y);
  vector[n] approx_y;
  vector[n] approx_var_y;
  vector[n+1] mode;
  vector[3 * n + 1] approx_results; // y, var, scaling, loglik
  real loglik = mode_[n+1];
  real diff = 1.0;
  int i = 0;
  //mode[1:n] = mode_[1:n]; 
  
  mode[1:n] = mode_[1:n] - xbeta; //adjust initial estimate with xbeta
  // check for bounds
  if (min(diagonal(Rt)) < 0.0) {
    reject("Negative standard deviation. ");
  }
  if (distribution != 2 && max(xbeta + mode[1:n]) > 50) {
    reject("Mean of the Poisson/negbin distribution > exp(50). ")
  }

  while(i < max_iter && diff > conv_tol) {

    vector[n+1] mode_new;
    approx_var_y = 1.0 ./ exp(xbeta + mode[1:n]);
      // note no xbeta here as it would be substracted in the smoother anyway
    approx_y = y .* approx_var_y + mode[1:n] - 1.0;
   
    mode_new = gaussian_smoother(approx_y, a1, P1, approx_var_y, Zt, Tt, Rt);
    // Problem with the approximation, potential divergence
    if (is_nan(mode_new[n+1]) || is_inf(mode_new[n+1]) ||
        (distribution != 2 && max(xbeta + mode_new[1:n]) > 50)) {
      reject("Error at iteration ", i, " of the approximation.");
    }
    diff = mean(square(mode_new[1:n] - mode[1:n]));
    mode = mode_new;
    loglik = mode[n+1];
    i = i + 1;
  }
  if (i == max_iter) {
    reject("Maximum number of iterations for approximation used. ");
  }

  approx_results[1:n] = approx_y;
  approx_results[(n+1):(2*n)] = approx_var_y;

  //mode term
  if(distribution == 1) {
    for(t in 1:n) {
      approx_results[2 * n + t] = y[t] * (xbeta[t] + mode[t]) - exp(xbeta[t] + mode[t]) +
        0.5 * ((approx_y[t] - mode[t])^2 / approx_var_y[t]);
    }
    approx_results[3 * n + 1] = loglik + sum(approx_results[(2*n+1): (3*n)]) +
       0.5 * sum(log(approx_var_y));
  } else {
    if (distribution == 2) {

    } else {

    }
  }
  
  return approx_results;
}

