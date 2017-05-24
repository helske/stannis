library(bssm)
data("exchange")

y <- as.numeric(dget("../mcmc_is/new_experiments/sp500.txt"))
model <- svm(y, rho = uniform(0.98,-0.9999,0.9999), 
  sd_ar = halfnormal(0.15, 5), mu = normal(0, 0, 10))

mcmc_out0 <- run_mcmc(model, n_iter = 5e4, nsim = 10, method = "pm")
mcmc_out <- run_mcmc(model, n_iter = 5e4, nsim = 10, method = "isc")
stan_inits <- list(list(mu = 0, phi = 0.98, sigma = 0.15))
stan_data <- list(y=y, n = length(y), initial_mode = c(model$initial_mode, -1e300))


fit <- stan("../stannis/exec/svm.stan", data = stan_data,
  init = stan_inits, chains = 1, iter = 10000)
fita <- stan("../stannis/exec/svm_approx.stan", data = stan_data,
  init = stan_inits, chains = 1, iter = 10000)


#####3
y <- as.numeric(dget("../mcmc_is/new_experiments/sp500.txt"))
model <- svm(y, rho = uniform(0.98,-0.9999,0.9999), 
  sd_ar = halfnormal(0.15, 5), mu = normal(0, 0, 10))
mcmc_out0 <- run_mcmc(model, n_iter = 5e4, nsim = 10, method = "pm")
mcmc_out <- run_mcmc(model, n_iter = 5e4, nsim = 10, method = "isc")
stan_inits <- list(list(mu = 0, phi = 0.98, sigma = 0.15))
stan_data <- list(y=y, n = length(y), initial_mode = c(model$initial_mode, -1e300))

  
fit <- stan("../stannis/exec/svm.stan", data = stan_data,
         init = stan_inits, chains = 1, iter = 10000)
# 
#     Elapsed Time: 407.39 seconds (Warm-up)
#     352.68 seconds (Sampling)
#     760.069 seconds (Total)
    
fita <- stan("../stannis/exec/svm_approx.stan", data = stan_data,
         init = stan_inits, chains = 1, iter = 10000)
    # Elapsed Time: 1262.31 seconds (Warm-up)
    # 1409.68 seconds (Sampling)
    # 2671.98 seconds (Total)
    # 
# > mcmc_out0$time
# user   system  elapsed 
# 1343.960    0.000 1357.90
# > mcmc_out$time
# user   system  elapsed 
# 1134.436    0.000 1134.794 

mcmc_out0
print(fit, c("phi", "sigma", "mu"))
# Inference for Stan model: svm.
# 1 chains, each with iter=10000; warmup=5000; thin=1; 
# post-warmup draws per chain=5000, total post-warmup draws=5000.
# 
# mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# phi    0.98       0 0.00  0.98  0.98  0.98  0.99  0.99  1719    1
# sigma  0.17       0 0.01  0.14  0.16  0.17  0.18  0.20  1443    1
# mu    -0.19       0 0.15 -0.49 -0.28 -0.19 -0.09  0.10  4238    1
# 
# Samples were drawn using NUTS(diag_e) at Sat Apr  8 10:52:16 2017.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#   convergence, Rhat=1).
print(fita, c("phi", "sd_x", "mu"))
# Inference for Stan model: svm_approx.
# 1 chains, each with iter=10000; warmup=5000; thin=1; 
# post-warmup draws per chain=5000, total post-warmup draws=5000.
# 
# mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# phi   0.98       0 0.00  0.98  0.98  0.98  0.99  0.99  2773    1
# sd_x  0.17       0 0.01  0.14  0.16  0.17  0.17  0.19  2659    1
# mu   -0.19       0 0.15 -0.48 -0.28 -0.19 -0.09  0.11  3169    1
# 
# Samples were drawn using NUTS(diag_e) at Sat Apr  8 11:36:57 2017.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#   convergence, Rhat=1).
mcmc_out0b <- run_mcmc(model, n_iter = 1e5, nsim = 10, method = "pm", delayed =TRUE)
mcmc_outb <- run_mcmc(model, n_iter = 1e5, nsim = 10, method = "isc")
