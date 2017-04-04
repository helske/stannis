library(rstan)
library(bssm)
library(stannis)
library(diagis)
set.seed(1)
n <- 1000
x <- cumsum(rnorm(n, 0, 0.2))
y <- rpois(n, exp(x))
ts.plot(y)

stan_data <- list(n = n, y = y, x1 = 0, P1 = 10)

stan_inits <- list(
  list(sd_x = 0.1, x = rep(0, n)),
  list(sd_x = 0.5, x = rep(1, n)),
  list(sd_x = 1, x = rep(-1, n)))

fit <- stan(file = 'stan/poisson_local_level.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 2000,  chains = 3, init = stan_inits)

## Approx
stan_inits <- list(
  list(sd_x = 0.1),
  list(sd_x = 0.5),
  list(sd_x = 1))

fita <- stan(file = 'stan/poisson_local_level_approx.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 2000,  chains = 3, init = stan_inits)


fita2 <- stan(file = 'stan/poisson_local_level_approx2.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 2000,  chains = 3, init = stan_inits)


c_time <- proc.time()
stan_out <- extract(fita, c("lp__", "sd_x", "jacobian"))
sd_x <- as.numeric(stan_out$sd_x)
lp__ <- as.numeric(stan_out$lp__ + stan_out$jacobian)
model <- ng_bsm(y, sd_level = halfnormal(0.1, 10), dist = "poisson", P1 = 10)
prior <- dnorm(sd_x, 0, 10, log = TRUE)
correction <- is_correction(model, sd_x, lp__, prior, 250)
c_time <- proc.time() - c_time


c_time <- proc.time()
stan_out <- extract(fita2, c("lp__", "sd_x", "jacobian", "approx_results"))
sd_x <- as.numeric(stan_out$sd_x)
lp__ <- as.numeric(stan_out$lp__ + stan_out$jacobian)
model <- ng_bsm(y, sd_level = halfnormal(0.1, 10), dist = "poisson", P1 = 10)
prior <- dnorm(sd_x, 0, 10, log = TRUE)
set.seed(1)
correction <- is_correction_psi(model, stan_out$approx_results, sd_x, lp__, stan_out$jacobian, 10)
set.seed(1)
correction2 <- is_correction_psi(model, stan_out$approx_results, sd_x, 
  as.numeric(stan_out$lp__), rep(0, length(sd_x)), 10)

c_time <- proc.time() - c_time
diagis::weight_plot(correction$weights)

###
c_time <- proc.time()
stan_out <- extract(fita, c("lp__", "sd_x", "jacobian"), perm = FALSE)
sd_x <- stan_out[,1,2]
lp__ <- as.numeric(stan_out[,1,1] + stan_out[,1,3])
model <- ng_bsm(y, sd_level = halfnormal(0.1, 10), dist = "poisson", P1 = 10)
prior <- dnorm(sd_x, 0, 10, log = TRUE)
correction <- is_correction(model, sd_x, lp__, prior, 250)
c_time <- proc.time() - c_time
effectiveSize(mcmc(extract(fit,"sd_x",perm=FALSE)[,,1]))
effectiveSize(mcmc(sd_x))
effectiveSize(mcmc(sd_x * correction$weights * length(sd_x)))
summary(mcmc(extract(fit,"sd_x",perm=FALSE)[,,1]))$stat
summary(mcmc(extract(fita,"sd_x",perm=FALSE)[,,1]))$stat
summary(mcmc(sd_x * correction$weights * length(sd_x)))$stat
##
diagis::weight_plot(correction$weights)
ess(correction$weights)
sum(get_elapsed_time(fit))
sum(get_elapsed_time(fita)) + c_time[3]
print(fit, "sd_x", digits = 4)
print(fita, "sd_x", digits = 4)
weighted_mean(sd_x, correction$weights)
weighted_se(sd_x, correction$weights)
# coda::spectrum0.ar(correction$weights * sd_x) / length(sd_x)
ts.plot(cbind(colMeans(extract(fit, "x")[[1]]), 
  rowMeans(correction$alpha[,1,]),
  weighted_mean(t(correction$alpha[,1,]), correction$weights)), col = 1:3)

