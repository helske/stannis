library(rstan)
set.seed(1)
n <- 100
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
print(fit, pars = "sd_x")

## Approx
stan_data <- list(n = n, y = as.numeric(y), x1 = 0, P1 = 10)

stan_inits <- list(
  list(sd_x = 0.1),
  list(sd_x = 0.5),
  list(sd_x = 1))

fita <- stan(file = 'stan/poisson_local_level_approx.stan', 
  data = stan_data, seed = 1,
  refresh = 0, iter = 50000,  chains = 1, init = stan_inits[1])
print(fita, pars = "sd_x")

stan_out <- extract(fita, c("lp__", "sd_x", "jacobian"))
sd_x <- stan_out$sd_x
lp__ <- stan_out$lp__
jacobian <- stan_out$jacobian

model <- ng_bsm(y, sd_level = halfnormal(0.1, 10), dist = "poisson", P1 = 10)
prior <- dnorm(sd_x, 0, 10, log = TRUE)
system.time(correction <- is_correction(model, sd_x, lp__ + jacobian, prior, 250))
hist(correction$weights)
hist(correction$posterior)
diagis::weight_plot(correction$weights)


w <- numeric(length(sd_x))
for (i in 1:length(sd_x)) {
  model$R[1] <- sd_x[i]
  w[i] <- logLik(model, nsim = 150, method ="bsf")
}

jac <- extract(fita,"jacobian")[[1]]
w <- numeric(length(sd_x))
for (i in 1:length(sd_x)) {
  model$R[1] <- sd_x[i]
  app <- gaussian_approx(model, max_iter = 1000,conv_tol=1e-15)
  aa <- c(smoother(app)$alpha)
  w[i] <- logLik(app) + dnorm(sd_x[i], 0, 10, log= TRUE) -lp__[i] - jac[i] +
   sum( dpois(y,exp(aa),log=TRUE)- dnorm(app$y,aa,app$H,log=TRUE))
}

library(diagis)
weighted_mean(as.numeric(sd_x), w)
mean(extract(fit, "sd_x")[[1]])
mean(extract(fita, "sd_x")[[1]])

sum(get_elapsed_time(fit))
sum(get_elapsed_time(fita))

sd(replicate(1000,logLik(model, nsim = 150, method ="bsf", seed =sample(1e8,size=1))))
