#' @importFrom rstan sampling
#' @importFrom bssm ng_bsm
#' @importFrom dplyr select starts_with
#' @export
stannis <- function(formula, data, level, slope, seasonal,
  distribution = "poisson", x1, P1, beta, iter = 7500, thin = 5, nsim_states = 10, 
  stan_inits, n_threads = 1, simplify_output = TRUE, is_seed = sample.int(.Machine$integer.max, 1), ...) {
  
  distr <- pmatch(distribution, "poisson") #pmatch(distribution, c("poisson", "binomial", "negative binomial"))
  
  # build xreg
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
  if (!missing(beta) && !is.matrix(beta)) beta <- matrix(beta, nrow = 1)
  ## local level or local linear trend model
  
  if (period == 1) {
    
    ## initial values
    if (missing(stan_inits)) {
      if (!missing(slope)) {
        stan_inits <- list(
          list(theta = c(0.01, 0.0001)), list(theta = c(0.1, 0.1)),
          list(theta = c(0.0001, 0.01)), list(theta = c(0.001, 0.001)))
      } else {
        stan_inits <- list(
          list(theta = c(0.01)), list(theta = c(0.1)),
          list(theta = c(0.0001)), list(theta = c(0.001)))
      }
    }
    if (k == 0) {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = if(!missing(slope)) halfnormal(slope[1],slope[2]),
        a1 = x1, P1 = P1, distribution = distribution)
      
      stan_data <- list(y = y, n = length(y),
        x1 = x1, P1 = P1, sd_prior_means = c(level[1], if(!missing(slope)) slope[1]),
        sd_prior_sds = c(level[2], if(!missing(slope)) slope[2]),
        initial_mode = c(model$initial_mode, -1e300), distribution = distr)
      
      fit <- sampling(
        if(missing(slope)) stannis:::stanmodels$ll_approx else stannis:::stanmodels$llt_approx,
        data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, cores = n_threads, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = if(!missing(slope)) halfnormal(slope[1],slope[2]),
        a1 = x1, P1 = P1, distribution = distribution, xreg = xreg,
        beta = normal(0, beta[, 1], beta[, 2]))
      
      stan_data <- list(y = y, n = length(y), k = k,
        x1 = x1, P1 = P1, sd_prior_means = c(level[1], if(!missing(slope)) slope[1]),
        sd_prior_sds = c(level[2], if(!missing(slope)) slope[2]),
        beta_prior_means = beta[, 1], beta_prior_sds = beta[, 2],
        initial_mode = c(model$initial_mode, -1e300), xreg = xreg, distribution = distr)
      if (k == 1) {
        dim(stan_data$beta_prior_means) <- dim(stan_data$beta_prior_sds) <- 1
      }
      fit <- sampling(stannis:::stanmodels$x_llt_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = c("Rt", "xbeta"), include = FALSE,
        iter = iter,  thin = thin, cores = n_threads, ...)
    }
    
  } else {
    # BSM, note we always include slope
    if (missing(stan_inits)) {
      stan_inits <- list(
        list(theta = rep(0.001, 3)), list(theta = c(0.01, 0.001, 0.001)),
        list(theta = c(0.001, 0.001, 0.1)),list(theta = c(0.01, 0.01, 0.001)))
    }
    if(k == 0L) {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = halfnormal(slope[1], slope[2]),
        sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = x1, P1 = P1,
        distribution = distribution)
      
      stan_data <- list(y = y, n = length(y), period = period,
        x1 = x1, P1 = P1, sd_prior_means = c(level[1], slope[1], seasonal[1]),
        sd_prior_sds = c(level[2], slope[2], seasonal[2]),
        initial_mode = c(model$initial_mode, -1e300), distribution = distr)
      
      fit <- sampling(stannis:::stanmodels$bsm_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, cores = n_threads, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = halfnormal(slope[1], slope[2]),
        sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = x1, P1 = P1,
        distribution = distribution, xreg = xreg,
        beta = normal(0, beta[, 1], beta[, 2]))
      
      
      
      stan_data <- list(y = y, n = length(y), k = k, period = period,
        x1 = x1, P1 = P1, sd_prior_means = c(level[1], slope[1], seasonal[1]),
        sd_prior_sds = c(level[2], slope[2], seasonal[2]),
        initial_mode = c(model$initial_mode, -1e300), distribution = distr,
        xreg = xreg, 
        beta_prior_means = beta[, 1], beta_prior_sds = beta[, 2])
      if (k == 1) {
        dim(stan_data$beta_prior_means) <- dim(stan_data$beta_prior_sds) <- 1
      }
      fit <- sampling(stannis:::stanmodels$x_bsm_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = c("Rt", "xbeta"), include = FALSE,
        iter = iter, thin = thin, ...)
    }
  }
  
  
  if (k==0) {
   
    stan_out <- as.data.frame(apply( 
      extract(fit, c("theta", "lp__", "jacobian", "approx_results"), permuted = FALSE), 
      3, identity))
    
    c_time <- proc.time()
    correction <- is_correction(model, as.matrix(select(stan_out, starts_with("approx_results"))),
      as.matrix(select(stan_out, starts_with("theta"))),
      select(stan_out, lp__)[, 1], select(stan_out, jacobian)[, 1], 
      nsim_states, n_threads, is_seed)
    c_time <- proc.time() - c_time
  } else {
    
    stan_out <- as.data.frame(apply( 
      extract(fit, c("theta", "beta", "lp__", "jacobian", "approx_results"), permuted = FALSE), 
      3, identity))
    
    c_time <- proc.time()
    correction <- is_correction(model, as.matrix(select(stan_out, starts_with("approx_results"))),
      as.matrix(cbind(select(stan_out, starts_with("theta")), select(stan_out, starts_with("beta")))), 
        select(stan_out, lp__)[, 1], select(stan_out, jacobian)[, 1], 
      nsim_states, n_threads, is_seed)
    c_time <- proc.time() - c_time
  }
  if (simplify_output) {
    list(theta = as.matrix(select(stan_out, starts_with("theta"))), 
      beta = if(k > 0) as.matrix(select(stan_out, starts_with("beta"))),
      states = correction$alpha, weights = c(correction$weights),
      s_time = sum(get_elapsed_time(fit)), c_time = c_time[3])
    
  } else {
  list(stan_fit = fit, theta = as.matrix(select(stan_out, starts_with("theta"))), 
    beta = if(k > 0) as.matrix(select(stan_out, starts_with("beta"))),
    correction = correction, s_time = get_elapsed_time(fit), c_time = c_time)
  }
}

