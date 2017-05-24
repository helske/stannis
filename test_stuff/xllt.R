set.seed(1)
n <- 500
k <- 3
x1 <- rnorm(n, mean = 1, sd = 0.1)
x2 <- cumsum(rnorm(n, sd = 0.1))
x3 <- runif(n, 0, 1)
xreg <- cbind(x1, x2, x3) #matrix(rnorm(mean = 1, n * k, sd = 0.1), n, k)
slope <- cumsum(c(0, rnorm(n-1, sd = 0.001)))
y <- rpois(n, exp(rowSums(xreg) + cumsum(slope + c(0, rnorm(n - 1, sd = 0.1)))))
ts.plot(y)



library(bssm)
model <- ng_bsm(y, xreg = xreg, beta = normal(rep(1, k), 0, 2), P1 = diag(100, 2),
  sd_level = halfnormal(0.1, 1), sd_slope = halfnormal(0.001, 1), distribution = "poisson")


out_da <- run_mcmc(model, n_iter = 50000, nsim = 10, delayed = TRUE)
out_is <- run_mcmc(model, n_iter = 50000, nsim = 10, method = "isc", const_m = FALSE)

dat <- as.data.frame(cbind(y, xreg))
names(dat) <- c("y",paste0("x",1:3))
bsm_out <- bsm(y~.,dat, distribution = "poisson", iter = 10000, thin=1,
  stan_inits = list(list(theta = c(0.1, 0.001), beta = rep(0, k))),
  beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
  sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(100, 2))


fit <- stan("exec/x_llt_poisson.stan", iter = 10000, thin=1, chains = 1,
  init = list(list(theta = c(0.1, 0.001), beta = rep(0, k), 
    slope_std = rep(0, n), level_std = rep(0, n))), data = list(y = as.integer(y), n = n, k =k, 
      xreg = xreg, beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
      sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(100, 2)))
