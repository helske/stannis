library(rstan)


stan_data <- list(n = length(Nile), y = Nile, x1 = 0, P1 = 1e4)

stan_inits <- list(
  list(sd_y = 100, sd_x = 100))

fit <- stan(file = 'stan/gaussian_ssm.stan', 
  data = stan_data, seed = 1,
  refresh = 1000, iter = 100,  chains = 1, init = stan_inits)
