library("bssm")
library("coda")
library("rstan")
library("diagis")
library("stannis")
library("doParallel")
library("foreach")

n <- 500
set.seed(1)
z1 <- runif(500, 0, 1)
z2 <- rexp(500, 5)
z3 <- cos(0.25 * (1:500))

slope <- cumsum(c(0, rnorm(500 - 1, sd = 0.001)))
level <- cumsum(slope + c(0, rnorm(500 - 1, sd = 0.1)))

y <- rpois(500, exp(z1 + z2 + z3 + level))

## standardize the covariates, Stan likes this...
sd_z1 <- sd(z1)
sd_z2 <- sd(z2)
sd_z3 <- sd(z3)
z1 <- (z1 - mean(z1)) / sd_z1
z2 <- (z2 - mean(z2)) / sd_z2
z3 <- (z3 - mean(z3)) / sd_z3


stan_data <- list(n = n, y = y, k = 3L, x1 = c(0, 0), P1 = diag(c(100, 1)), 
  sd_prior_means = rep(0, 2), sd_prior_sds = rep(1, 2),
  xreg = cbind(z1, z2, z3), beta_prior_means = rep(0, 3), beta_prior_sds = rep(10, 3))

stan_inits <- list(list(theta = c(0.1, 0.001), beta = c(sd_z1, sd_z2, sd_z3), 
  level_std = rep(1, n), slope_std = rep(0, n)))

res_init_states_std <- sampling(stannis:::stanmodels$x_llt_poisson, 
    data = stan_data, refresh = 1000, seed = 1,
    iter = 5e4, chains = 1, cores = 1, init = stan_inits)

stan_inits <- list(list(theta = c(0.1, 0.001), beta = c(sd_z1, sd_z2, sd_z3)))

res_noinit_states_std <- sampling(stannis:::stanmodels$x_llt_poisson, 
  data = stan_data, refresh = 1000, seed = 1,
  iter = 5e4, chains = 1, cores = 1, init = stan_inits)


set.seed(1)
z1 <- runif(500, 0, 1)
z2 <- rexp(500, 5)
z3 <- cos(0.25 * (1:500))

slope <- cumsum(c(0, rnorm(500 - 1, sd = 0.001)))
level <- cumsum(slope + c(0, rnorm(500 - 1, sd = 0.1)))

y <- rpois(500, exp(z1 + z2 + z3 + level))

stan_data <- list(n = n, y = y, k = 3L, x1 = c(0, 0), P1 = diag(c(100, 1)), 
  sd_prior_means = rep(0, 2), sd_prior_sds = rep(1, 2),
  xreg = cbind(z1, z2, z3), beta_prior_means = rep(0, 3), beta_prior_sds = rep(10, 3))

stan_inits <- list(list(theta = c(0.1, 0.001), beta = c(1, 1, 1)))
res_noinit_states_nostd <- sampling(stannis:::stanmodels$x_llt_poisson, 
  data = stan_data, refresh = 1000, seed = 1,
  iter = 5e4, chains = 1, cores = 1, init = stan_inits)

stan_inits <- list(list(theta = c(0.1, 0.001), beta = c(1, 1, 1), std_level = rep(1, n), std_slope = rep(0, n)))
res_init_states_nostd <- sampling(stannis:::stanmodels$x_llt_poisson, 
  data = stan_data, refresh = 1000, seed = 1,
  iter = 5e4, chains = 1, cores = 1, init = stan_inits)


print(res_init_states_std, pars = c("theta", "beta"))
print(res_noinit_states_std, pars = c("theta", "beta"))
print(res_init_states_nostd, pars = c("theta", "beta"))
print(res_noinit_states_nostd, pars = c("theta", "beta"))
