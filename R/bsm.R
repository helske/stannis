#' @export
bsm <- function(formula, data, distribution, x1, P1, beta_prior_means, beta_prior_sds, 
  sd_prior_means, sd_prior_sds, iter = 7500, thin = 5, 
  nsim_states = 10, stan_inits, ...) {
  
  distribution <- pmatch(distribution, c("poisson", "binomial", "negative binomial"))
  
  if(inherits(formula, "formula")) {
    if (missing(data)) {
      data <- environment(formula)
      period <- NULL
    } else {
      period <- tsp(data)[3]
      data <- as.data.frame(data)
    }
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$na.action <- as.name("na.pass")
    mf <- eval(mf, parent.frame())
    # remove the intercept
    attr(attr(mf, "terms"), "intercept") <- 0L
    y <- model.response(mf, "numeric")
    xreg <- model.matrix(attr(mf, "terms"), mf)
    k <- ncol(xreg)
    if (is.null(period)) {
      period <- frequency(y)
    }
  } else {
    y <- formula
    k <- 0L
    period <- frequency(y)
  }
  
  if (period == 1) {
    
    if (missing(stan_inits)) {
      if (k == 0) {
        stan_inits <- list(
          list(theta = c(0.01, 0.0001)), list(theta = c(0.1, 0.1)),
          list(theta = c(0.0001, 0.01)),list(theta = c(0.001, 0.001)))
      } else {
        stan_inits <- list(
          list(theta = c(0.01, 0.0001), beta = rep(0, k)), 
          list(theta = c(0.1, 0.1), beta = rep(0, k)),
          list(theta = c(0.0001, 0.01), beta = rep(0, k)),
          list(theta = c(0.001, 0.001), beta = rep(0, k)))
      }
    }
    if (k == 0) {
      model <- ng_bsm(y, sd_level = halfnormal(sd_prior_means[1], sd_prior_sds[1]),
        sd_slope = halfnormal(sd_prior_means[2], sd_prior_sds[2]),
        a1 = x1, P1 = P1, distribution = distribution)
      
      stan_data <- list(y = y, n = length(y), period = period,
        x1 = x1, P1 = P1, sd_prior_means = sd_prior_means, sd_prior_sds = sd_prior_sds,
        initial_mode = c(model$initial_mode, -1e300))
      
      fit <- sampling(stannis:::stanmodels$llt_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(sd_prior_means[1], sd_prior_sds[1]),
        sd_slope = halfnormal(sd_prior_means[2], sd_prior_sds[2]),
        a1 = x1, P1 = P1, distribution = distribution, xreg = xreg, 
        beta = normal(0, beta_prior_means, beta_prior_sds))
      
      stan_data <- list(y = y, n = length(y), k = k, period = period,
        x1 = x1, P1 = P1, sd_prior_means = sd_prior_means, sd_prior_sds = sd_prior_sds,
        beta_prior_means = beta_prior_means, beta_prior_sds = beta_prior_sds,
        initial_mode = c(model$initial_mode, -1e300), xreg = xreg)
      
      fit <- sampling(stannis:::stanmodels$x_llt_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = c("Rt", "xbeta"), include = FALSE,
        iter = iter,  thin = thin, ...)
    }
    
  } else {
    
    if (missing(stan_inits)) {
      if (k == 0L) {
        stan_inits <- list(
          list(theta = rep(0.001, 3)), list(theta = c(0.01, 0.001, 0.001)),
          list(theta = c(0.001, 0.001, 0.1)),list(theta = c(0.01, 0.01, 0.001)))
      } else {
        stan_inits <- list(
          list(theta = rep(0.001, 3), beta = rep(0, k)),
          list(theta = c(0.01, 0.001, 0.001), beta = rep(0, k)),
          list(theta = c(0.001, 0.001, 0.1), beta = rep(0, k)),
          list(theta = c(0.01, 0.01, 0.001), beta = rep(0, k)))
      }
    }
    if(k == 0L) {
      model <- ng_bsm(y, sd_level = halfnormal(sd_prior_means[1], sd_prior_sds[1]),
        sd_slope = halfnormal(sd_prior_means[2], sd_prior_sds[2]),
        sd_seasonal = halfnormal(sd_prior_means[3], sd_prior_sds[3]),  a1 = x1, P1 = P1,
        distribution = distribution)
      
      stan_data <- list(y = y, n = length(y), period = period,
        x1 = x1, P1 = P1, sd_prior_means = sd_prior_means, sd_prior_sds = sd_prior_sds,
        initial_mode = c(model$initial_mode, -1e300))
      
      fit <- sampling(stannis:::stanmodels$bsm_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(sd_prior_means[1], sd_prior_sds[1]),
        sd_slope = halfnormal(sd_prior_means[2], sd_prior_sds[2]),
        sd_seasonal = halfnormal(sd_prior_means[3], sd_prior_sds[3]),  a1 = x1, P1 = P1,
        distribution = distribution, xreg = xreg, 
        beta = normal(0, beta_prior_means, beta_prior_sds))
      
      stan_data <- list(y = y, n = length(y), k = k, period = period,
        x1 = x1, P1 = P1, sd_prior_means = sd_prior_means, sd_prior_sds = sd_prior_sds,
        beta_prior_means = beta_prior_means, beta_prior_sds = beta_prior_sds,
        initial_mode = c(model$initial_mode, -1e300), xreg = xreg)
      
      fit <- sampling(stannis:::stanmodels$x_bsm_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits, 
        pars = c("Rt", "xbeta"), include = FALSE,
        iter = iter, thin = thin, ...)
    }
  }
  
  
  if (k==0) {
    c_time <- proc.time()
    stan_out <- extract(fit, c("lp__", "theta", "jacobian", "approx_results"))
    correction <- is_correction_psi(model, stan_out$approx_results,
      stan_out$theta, stan_out$lp__, stan_out$jacobian, nsim_states)
    c_time <- proc.time() - c_time
  } else {
    c_time <- proc.time()
    stan_out <- extract(fit, c("lp__", "theta", "beta", "jacobian", "approx_results"))
    correction <- is_correction_psi(model, stan_out$approx_results,
      cbind(stan_out$theta, stan_out$beta), stan_out$lp__, stan_out$jacobian, nsim_states)
    c_time <- proc.time() - c_time
  }
  list(stan_fit = fit, correction = correction, s_time = get_elapsed_time(fit), c_time = c_time)
}

