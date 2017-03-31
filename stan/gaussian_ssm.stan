functions {
 
 real local_level_lpdf(vector y, real x1, real P1, real sd_y, real sd_x) {
   
   real LOG2PI = log(2 * pi());
   
   int n = rows(y);
   real loglik = 0.0;
   real x = x1;
   real P = P1;
   real var_y = sd_y * sd_y;
   real var_x = sd_x * sd_x;
   
   for(t in 1:n) {
    real v = y[t] - x;
    real F = P + var_y;
    real K = P / F;
    x = x + K * v;
    P = P * (1.0 - K) + var_x;
    loglik = loglik - 0.5 * (LOG2PI + log(F) + v * v / F);
   }  
   return loglik;
 }
 
}
data {
  int<lower=0> n;
  vector[n] y;
  real x1;
  real P1;
}

parameters {
  real<lower=0> sd_y;
  real<lower=0> sd_x;
}

model {
  target += local_level_lpdf(y | x1, P1, sd_y, sd_x);
}
