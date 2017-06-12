library("bssm")
library("coda")
library("rstan")
library("diagis")
library("stannis")
library("doParallel")
library("foreach")

ire_experiment_llt <- function(n_iter,
  nsim_states = 10, seed = sample(.Machine$integer.max, 1), method){
  
  # simulate the data (with fixed seed)
  set.seed(123)
  n <- 100
  
  slope <- cumsum(c(0, rnorm(n - 1, sd = 0.01)))
  level <- cumsum(slope + c(0, rnorm(n - 1, sd = 0.01)))
  y <- rpois(n, exp(level))
  
  
  results <- data.frame(method = method, 
    time = 0, "theta_1" = 0, "theta_2" = 0, 
    "level_1" = 0, "slope_1" = 0, 
    "level_n" = 0, "slope_n" = 0, 
    "divergent" = 0, "treedepth" = 0)
  
  set.seed(seed)
  
  if(method == "stannis") {
    
    res <- stannis(y, iter = n_iter, 
      level = c(0, 1), slope = c(0, 0.1), refresh = 0, a1 = c(0, 0), P1 = diag(c(10, 0.1)),
      stan_inits = list(list(theta = c(0.01, 0.01))))
    
    results[1, 2] <- res$stan_time + res$correction_time
    results[1, 3:4] <- res$mean_theta
    results[1, 5:6] <- weighted_mean(t(res$states[1,,]), res$weights)
    results[1, 7:8] <- weighted_mean(t(res$states[n,,]), res$weights)
    saveRDS(results, file = paste0("stannis_llt_preresults100_",seed,".rda"))
    
    return(results)
  }
  if(method == "Stan") {
    stan_data <- list(n = n, y = y, a1 = c(0, 0), P1 = diag(c(10, 0.1)), 
      sd_prior_means = rep(0, 2), sd_prior_sds = c(1, 0.1))
    stan_inits <- list(list(theta = c(0.01, 0.01), 
      slope_std = rep(0, n), level = rep(0, n)))
    
    res <- sampling(stannis:::stanmodels$llt_poisson, 
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      data = stan_data, refresh = 0, 
      iter = n_iter, chains = 1, cores = 1, init = stan_inits)
    results[1, 2] <- sum(get_elapsed_time(res))
    results[1, 3:4] <- summary(res, pars = "theta")$summary[, "mean"]
    results[1, 5] <- summary(res, pars = "level[1]")$summary[, "mean"]
    results[1, 6] <- summary(res, pars = "slope[1]")$summary[, "mean"]
    results[1, 7] <- summary(res, pars = paste0("level[",n,"]"))$summary[, "mean"]
    results[1, 8] <- summary(res, pars = paste0("slope[",n,"]"))$summary[, "mean"]
    diags <- get_sampler_params(res, inc_warmup = FALSE)[[1]]
    results[1, 9] <- sum(diags[, "divergent__"])
    results[1, 10] <- sum(diags[, "treedepth__"] >= 15)
    saveRDS(results, file = paste0("stan_llt_preresults100_",seed,".rda"))
    
    return(results)
  }
  
  model <- ng_bsm(y, P1 = diag(c(10, 0.1)),
    sd_level = halfnormal(0.01, 1), 
    sd_slope = halfnormal(0.01, 0.1), distribution = "poisson")
  if(method == "isc") {
    res <- run_mcmc(model, n_iter = n_iter, nsim_states = 10, 
      method = "isc", n_threads = 1) 
    results[1, 2] <- res$time[3]
    results[1, 3:4] <- weighted_mean(res$theta, res$weights * res$counts)
    results[1, 5:6] <- weighted_mean(t(res$alpha[1,,]), res$weights * res$counts)
    results[1, 7:8] <- weighted_mean(t(res$alpha[n,,]), res$weights * res$counts)
    
    return(results)
  }
  if(method == "da") {
    res <- run_mcmc(model, n_iter = n_iter, nsim_states = 10, method = "pm") 
    results[1, 2] <- res$time[3]
    results[1, 3:4] <- weighted_mean(res$theta, res$counts)
    results[1, 5:6] <- weighted_mean(t(res$alpha[1,,]), res$counts)
    results[1, 7:8] <- weighted_mean(t(res$alpha[n,,]), res$counts)
    
    return(results)
  }
  
  res <- run_mcmc(model, n_iter = n_iter, nsim_states = 10, method = "pm", 
    delayed_acceptance = FALSE) 
  results[1, 2] <- res$time[3]
  results[1, 3:4] <- weighted_mean(res$theta, res$counts)
  results[1, 5:6] <- weighted_mean(t(res$alpha[1,,]), res$counts)
  results[1, 7:8] <- weighted_mean(t(res$alpha[n,,]), res$counts)
  results
}


cl<-makeCluster(16)
registerDoParallel(cl)

results <- 
  foreach (i = 1:100, .combine=rbind, .packages = c("bssm", "diagis", "stannis", "rstan")) %dopar% 
  ire_experiment_llt(n_iter = 4e4, seed = i, method = "Stan")
saveRDS(results, file = "stan_llt_iter4e4.rda")


results <- 
  foreach (i = 1:100, .combine=rbind, .packages = c("bssm", "diagis", "stannis", "rstan")) %dopar% 
  ire_experiment_llt(n_iter = 4e4, seed = i, method = "stannis")
saveRDS(results, file = "stannis_llt_iter4e4.rda")


results <- 
  foreach (i = 1:100, .combine=rbind, .packages = c("bssm", "diagis", "stannis", "rstan")) %dopar% 
  ire_experiment_llt(n_iter = 4e4, seed = i, method = "isc")
saveRDS(results, file = "is_llt_iter4e4.rda")

results <- 
  foreach (i = 1:100, .combine=rbind, .packages = c("bssm", "diagis", "stannis", "rstan")) %dopar% 
  ire_experiment_llt(n_iter = 4e4, seed = i, method = "da")
saveRDS(results, file = "da_llt_iter4e4.rda")

results <- 
  foreach (i = 1:100, .combine=rbind, .packages = c("bssm", "diagis", "stannis", "rstan")) %dopar% 
  ire_experiment_llt(n_iter = 4e4, seed = i, method = "pm")
saveRDS(results, file = "pm_llt_iter4e4.rda")



stopCluster(cl)

## reference values

set.seed(123)
n <- 100

slope <- cumsum(c(0, rnorm(n - 1, sd = 0.01)))
level <- cumsum(slope + c(0, rnorm(n - 1, sd = 0.01)))
y <- rpois(n, exp(level))

model <- ng_bsm(y, P1 = diag(c(10, 0.1)),
  sd_level = halfnormal(0.01, 1), 
  sd_slope = halfnormal(0.01, 0.1), distribution = "poisson")
res <- run_mcmc(model, n_iter = 1.1e6, n_burnin = 1e5, delayed_acceptance = FALSE, nsim = 10)

theta <- weighted_mean(res$theta, res$counts)
level_1 <- weighted_mean(res$alpha[1,1,], res$counts)
level_n <- weighted_mean(res$alpha[n,1,], res$counts)
slope_1 <- weighted_mean(res$alpha[1,2,], res$counts)
slope_n <- weighted_mean(res$alpha[n,2,], res$counts)
reference_llt <- c(theta_1 = theta[1], theta_2 = theta[2],
  level_1 = level_1, level_n = level_n, 
  slope_1 = slope_1, slope_n = slope_n)
saveRDS(reference_llt, file = "reference_llt.rda")
