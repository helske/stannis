library(rstan)
set.seed(123)
n <- 1000
x <- cumsum(rnorm(n, 0, 0.1))
y <- rpois(n, exp(x))
ts.plot(y)

stan_data <- list(n = n, y = y, x1 = 0, P1 = 10)

stan_inits <- list(
  list(sd_x = 0.1, x = rep(0, n)),
  list(sd_x = 0.5, x = rep(1, n)),
  list(sd_x = 1, x = rep(-1, n)))

fit <- stan(file = 'stan/poisson_local_level.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 50000,  chains = 3, init = stan_inits)

## Approx
stan_inits <- list(
  list(sd_x = 0.1),
  list(sd_x = 0.5),
  list(sd_x = 1))

fita <- stan(file = 'stan/poisson_local_level_approx.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 50000,  chains = 3, init = stan_inits)

c_time <- proc.time()
stan_out <- extract(fita, c("lp__", "sd_x", "jacobian"))
sd_x <- stan_out$sd_x
lp__ <- stan_out$lp__ + stan_out$jacobian
model <- ng_bsm(y, sd_level = halfnormal(0.1, 10), dist = "poisson", P1 = 10)
prior <- dnorm(sd_x, 0, 10, log = TRUE)
correction <- is_correction(model, sd_x, lp__, prior, 250)
c_time <- proc.time() - c_time

diagis::weight_plot(correction$weights)
sum(get_elapsed_time(fit))
sum(get_elapsed_time(fita)) + c_time[3] / 1000
print(fit, "sd_x", digits = 4)
print(fita, "sd_x", digits = 4)
weighted_mean(sd_x, correction$weights)
weighted_se(sd_x, correction$weights)
# coda::spectrum0.ar(correction$weights * sd_x) / length(sd_x)
ts.plot(cbind(rowMeans(extract(fit, "x")), 
  weighted_mean(correction$alpha[,1,], correction$weights)), col = 1:2)

