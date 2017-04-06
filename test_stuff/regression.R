set.seed(1)
n <- 100
k <- 3
x1 <- rnorm(n, mean = 1, sd = 0.1)
x2 <- cumsum(rnorm(n, sd = 0.1))
x3 <- runif(n, 0, 1)
xreg <- cbind(x1, x2, x3) #matrix(rnorm(mean = 1, n * k, sd = 0.1), n, k)
slope <- cumsum(c(0, rnorm(n-1, sd = 0.001)))
y <- rpois(n, exp(rowSums(xreg) + cumsum(slope + c(0, rnorm(n - 1, sd = 0.1)))))
ts.plot(y)

library(bssm)
model <- ng_bsm(y, xreg = xreg, beta = normal(rep(0, k), 0, 2), P1 = diag(100, 2),
  sd_level = halfnormal(0.1, 1), sd_slope = halfnormal(0.001, 1), distribution = "poisson")

out <- run_mcmc(model, n_iter = 1e4, nsim = 10, method = "isc")
out2 <- run_mcmc(model, n_iter = 1e4, nsim = 10, method = "pm")

dat <- as.data.frame(cbind(y, xreg))
names(dat) <- c("y",paste0("x",1:3))
bsm_out <- stannis::bsm(y~.,dat, distribution = "poisson", iter = 2000, thin=1,
  stan_inits = list(list(theta = c(0.1, 0.001), beta = rep(0, k))),
  beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
  sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(100, 2))
print(bsm_out$stan_fit, "beta")

summary(mcmc(extract(bsm_out$stan_fit, "beta",perm=FALSE)[,1,1]))

# Iterations = 1:5000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 5000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 2.057807       0.563159       0.007964       0.007575 
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5% 
# 0.9591 1.6724 2.0467 2.4312 3.1834 
out <- run_mcmc(model, n_iter = 5000, nsim = 10)

dat <- as.data.frame(cbind(y, xreg))
names(dat) <- c("y",paste0("x",1:3))
bsm_out <- stannis::bsm(y~.,dat, distribution = "poisson", iter = 2000, thin=1,
  stan_inits = list(list(theta = c(0.1, 0.001), beta = rep(0, k))),
  beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
  sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(100, 2))
