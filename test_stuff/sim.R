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
model <- ng_bsm(y, xreg = xreg, beta = normal(rep(1, k), 0, 2), P1 = diag(100, 2),
  sd_level = halfnormal(0.1, 1), sd_slope = halfnormal(0.001, 1), distribution = "poisson")

dat <- as.data.frame(cbind(y, xreg))
names(dat) <- c("y",paste0("x",1:3))

nsim <- 100
mean_is <- mean_stan <- matrix(NA, 5, nsim)
time_is <- time_stan <- numeric(nsim)
for(i in 1:nsim){
  
  bsm_out <- bsm(y~.,dat, distribution = "poisson", iter = 2000, thin=1, refresh = 0,
    stan_inits = list(list(theta = c(0.1, 0.001), beta = rep(0, k))),
    beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
    sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(c(25, 1), 2))
  
  output <- extract(bsm_out$stan_fit, c("beta", "theta"))
  mean_is[,i] <- c(weighted_mean(output$theta, bsm_out$correction$weights), 
    weighted_mean(output$beta, bsm_out$correction$weights))
  time_is[i] <- sum(bsm_out$s_time) + bsm_out$c_time[3]
  fit <- stan("exec/x_llt_poisson.stan", iter = 2000, thin=1, chains = 1, refresh = 0,
    init = list(list(theta = c(0.1, 0.001), beta = rep(0, k), 
      slope_std = rep(0, n), level_std = rep(0, n))), data = list(y = as.integer(y), n = n, k =k, 
        xreg = xreg, beta_prior_means = rep(0, k), beta_prior_sds = rep(2, k), 
        sd_prior_means = c(0, 0), sd_prior_sds = c(1, 1), x1 = c(0,0), P1 = diag(c(25, 1), 2)))
  output <- extract(fit, c("beta", "theta"))
  mean_stan[,i] <- c(colMeans(output$theta), colMeans(output$beta))
  time_stan[i] <- sum(get_elapsed_time(fit))
  apply(mean_is, 1, sd, na.rm = TRUE)
  apply(mean_stan, 1, sd, na.rm = TRUE)
  print(i)
}
