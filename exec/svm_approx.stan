functions {

vector svm_smoother(vector y, real mu, real phi, real sd_x, vector var_y) {

  int n = rows(y);
  real loglik = 0.0;
  real x = mu;
  real P = sd_x * sd_x / sqrt(1 - phi * phi) ;
  real var_x = sd_x * sd_x;

  vector[n] v;
  vector[n] F;
  vector[n] K;
  vector[n+1] r;

  for (t in 1:n) {
    F[t] = P + var_y[t];
    if (F[t] > 1.0e-8) {
      v[t] = y[t] - x;
      K[t] = P / F[t];
      x = phi * (x + K[t] * v[t]);
      P = phi * phi * P * (1.0 - K[t]) + var_x;
      loglik = loglik - 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = phi * x;
      P = phi * phi * P + var_x;
    }
  }  
  
  r[n+1] = 0.0;
  for (tt in 1:n) {
    int t = n + 1 - tt;
    if(F[t] > 1.0e-8) {
      r[t] = v[t] / F[t] + phi * (1.0 - K[t]) * r[t+1];
    } else {
      r[t] = phi * r[t+1];
    }
  }

  r[1] = mu + sd_x * sd_x / sqrt(1 - phi * phi) * r[1];
  for (t in 2:n) {
    r[t] = phi * r[t-1] + var_x * r[t];
  }
  r[n+1] = loglik;
  return r;
  
}

 vector svm_approx(vector y, vector iy_nz2, real mu, real phi, real sd_x, vector mode_) {

  int n = rows(y);
  vector[n] approx_y;
  vector[n] approx_var_y;
  vector[n+1] mode = mode_;
  vector[n+1] mode_new;
  vector[3 * n + 1] approx_results; // y, var, scaling, loglik
  real loglik = mode[n+1];
  real diff = 1.0;
  int i = 1;
   
  while(i < 25 && diff > 1.0e-8) {
    
    approx_var_y = 2.0 * exp(mode[1:n]) .* iy_nz2;
    approx_y = mode[1:n] + 1.0 - 0.5 * approx_var_y;
  
    mode_new = svm_smoother(approx_y, mu, phi, sd_x, approx_var_y);
    // Problem with the approximation, potential divergence
    if (is_nan(mode_new[n+1]) || is_inf(mode[n+1])) {
      reject("Error at iteration ", i, " of the approximation.");
    }
    diff = mean(square(mode_new[1:n] - mode[1:n]));
    mode = mode_new;
    loglik = mode[n+1];
    i = i + 1;
  }
  if (i == 25) {
    reject("Maximum number of iterations for approximation used. ");
  }
  
  approx_results[1:n] = approx_y;
  approx_results[(n+1):(2*n)] = approx_var_y;
  
  //mode term
 
  for(t in 1:n) {
    approx_results[2 * n + t] = -0.5 * (mode[t] + y[t]^2 * exp(-mode[t])) +
      0.5 * (approx_y[t] - mode[t])^2 / approx_var_y[t] + log(approx_var_y[t]);
  }

  approx_results[3 * n + 1] = loglik + sum(approx_results[(2*n+1): (3*n)]);
  return approx_results;
}

}

data {
  int<lower=0> n;
  vector[n] y;
  vector[n+1] initial_mode;
}

transformed data {
  vector[n] iy_nz2 = y;
  for (t in 1:n) {
    if (fabs(y[t]) < 1.0e-4) {
      iy_nz2[t] = 1.0e-4;
    }
    iy_nz2[t] = 1.0 / iy_nz2[t]^2;
  }
}

parameters {
  real mu;
  real<lower=-1,upper=1> phi; 
  real<lower=0> sd_x;
}


transformed parameters {
  vector[3 * n + 1] approx_results;
  approx_results =  svm_approx(y, iy_nz2, mu, phi, sd_x, initial_mode);
}

model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(sd_x | 0, 5);
  target += approx_results[3 * n + 1];
}
