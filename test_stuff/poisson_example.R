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
model <- ng_bsm(y, xreg = xreg, beta = normal(rep(1, k), 0, 2), P1 = diag(1, 2),
  sd_level = halfnormal(0.1, 1), sd_slope = halfnormal(0.001, 1), distribution = "poisson")

sd(replicate(100, logLik(model, nsim = 10, method="psi", seed = sample(1e8, size=1))))

out_da <- run_mcmc(model, n_iter = 5000, nsim = 10, delayed = TRUE)
out_is <- run_mcmc(model, n_iter = 5000, nsim = 10, method = "isc")
out_bsf <- run_mcmc(model, n_iter = 5000, nsim = 1000, method = "pm",delayed=TRUE,simulation="bsf")
dat <- as.data.frame(cbind(y, xreg))
names(dat) <- c("y",paste0("x",1:3))
bsm_out <- bsm(y~.,dat, distribution = "poisson", iter = 2000, thin=1,
  stan_inits = list(list(theta = c(0.1, 0.001), beta = rep(0, k))),
  beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
  sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(1, 2))
