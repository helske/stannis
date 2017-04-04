functions {
  
 vector local_level_smoother(vector y, real x1, real P1, vector var_y, real sd_x) {
   
   real LOG2PI = log(2 * pi());
   
   int n = rows(y);
   real loglik = 0.0;
   real x = x1;
   real P = P1;
   real var_x = sd_x * sd_x;
   
   vector[n] v;
   vector[n] F;
   vector[n] K;
   vector[n+1] r;
   
   for (t in 1:n) {
    v[t] = y[t] - x;
    F[t] = P + var_y[t];
    K[t] = P / F[t];
    x = x + K[t] * v[t];
    P = P * (1.0 - K[t]) + var_x;
    loglik = loglik - 0.5 * (LOG2PI + log(F[t]) + v[t] * v[t] / F[t]);
   }  
   
   
   r[n+1] = 0.0;
   for (tt in 1:n) {
     int t = n + 1 - tt;
     r[t] = v[t] / F[t] + (1.0 - K[t]) * r[t+1];
   }
   
   r[1] = x1 + P1 * r[1];
   for (t in 2:n) {
     r[t] = r[t-1] + var_x * r[t];
   }
   r[n+1] = loglik;
   return r;
 }
 
 vector poisson_local_level_approx(int[] y, real x1, real P1, real sd_x, vector theta_) {

   int n = size(y);
   vector[n] yreal;
   vector[n] approx_y;
   vector[n] approx_var_y;
   vector[n+1] theta = theta_;
   vector[3 * n + 1] approx_results; // y, var, scaling, loglik
   real loglik = theta[n+1];
   real diff = 1;
   int i=1;
   for(t in 1:n) {
     yreal[t] = y[t];
   }
   while(i < 100 && diff > 1.0e-8) {
     
    approx_var_y = 1.0 ./ exp(theta[1:n]);
    approx_y = yreal .* approx_var_y + theta[1:n] - 1.0;
    theta = local_level_smoother(approx_y, x1, P1, approx_var_y, sd_x);
    diff = fabs(theta[n+1] - loglik) / (0.1 + fabs(theta[n+1]));
    loglik = theta[n+1];
    if (is_nan(loglik) || is_inf(loglik)) {
      // print("sd_x: ", sd_x);
      // print("approx_var_y: ", approx_var_y);
      // print("approx_y: ", approx_y);
      // print("theta: ", theta);
      reject("Approximation failed due to non-finite loglikelihood at iteration ", i);
    }
    i = i + 1;
   }
   if (i == 100) {
     reject("Maximum number of iterations for approximation used. ");
   }
  
   approx_results[1:n] = approx_y;
   approx_results[(n+1):(2*n)] = approx_var_y;
   //mode term
   for(t in 1:n) {
      approx_results[2 * n + t] = y[t] * theta[t] - exp(theta[t]) +
     0.5 * (approx_y[t] - theta[t])^2 / approx_var_y[t] + log(approx_var_y[t]);
   }

   approx_results[3 * n + 1] = loglik + sum(approx_results[(2*n+1): (3*n)]);
   return approx_results;
 }
 
}

data {
  int<lower=0> n;
  int<lower=0> y[n];
  real x1;
  real P1;
}
transformed data {

  vector[n+1] theta_init;
  theta_init[n+1] = -1e300; //store loglik here
  for(t in 1:n) {
    if(y[t]<0.1) {
      theta_init[t] = log(0.1);
    } else {
      theta_init[t] = log(y[t]);
    }
  }

}
parameters {
  real<lower=0> sd_x;
}
transformed parameters {
    vector[3 * n + 1] approx_results = poisson_local_level_approx(y, x1, P1, sd_x, theta_init);
}
model {
  target += normal_lpdf(sd_x | 0, 10);
  target += approx_results[3 * n + 1];
}

generated quantities {
  real jacobian;
  jacobian = -log(sd_x);
}
