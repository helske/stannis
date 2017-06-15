#' Structural Poisson Time Series with External Covariates
#' 
#' \code{stannis} performs fully Bayesian inference of structural Poisson time series optionally with covariates. 
#' 
#' \code{stannis} first obtains approximate marginal posterior using Markov chain Monte Carlo via 
#' Hamiltonian Monte Carlo sampler of \code{Stan}, and then performs importance sampling type correction of this 
#' posterior sample in order to obtain "exact" joint posterior. The correction step is done using particle filtering.
#'
#' @import methods
#' @importFrom rstan sampling get_elapsed_time extract
#' @importFrom stats model.matrix model.response tsp frequency
#' @importFrom Rcpp loadModule
#' @importFrom bssm ng_bsm halfnormal normal smoother
#' @importFrom dplyr select_ starts_with
#' @importFrom coda spectrum0.ar
#' @importFrom diagis weighted_mean weighted_var
#' @param formula Object of class \code{formula} which defines the regression part. 
#' Or a single vector of observations in case of no covariates.
#' @param data Optional data frame (or object coercible to such) containing the data.
#' @param level A vector of length two, defining the mean and standard deviation of the 
#' truncated Gaussian prior for the standard deviation of the level noise term.
#' @param slope A vector of length two, defining the mean and standard deviation of the 
#' truncated Gaussian prior for the standard deviation of the slope noise term.
#' @param seasonal A vector of length two, defining the mean and standard deviation of the 
#' truncated Gaussian prior for the standard deviation of the seasonal noise term.
#' @param a1 Prior mean for the state vector at time 1.
#' @param P1 Prior covariance matrix for state vector at time 1.
#' @param beta A matrix with \eqn{k} rows and 2 columns, where first columns defines the 
#' prior means of the Gaussian priors of the corresponding \eqn{k} regression coefficients, 
#' and the second column defines the the standard deviations of those prior distributions.
#' @param iter Number of iterations for the MCMC.
#' @param thin Thinning interval used in MCMC.
#' @param nsim_states Number of particles used in particle filter.
#' @param stan_inits Initial values for the MCMC as a list of lists.
#' @param n_threads Number of threads to use in parallel.
#' @param simplify_output If \code{TRUE} (default), does not return the original output from \code{Stan}.
#' @param is_seed Seed for the importance sampling correction.
#' @export
stannis <- function(formula, data, level, slope, seasonal,
 a1, P1, beta, iter, thin = 1, nsim_states = 10, 
  stan_inits, n_threads = 1, simplify_output = TRUE, is_seed = sample.int(.Machine$integer.max, 1), ...) {
  
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
          list(theta = c(0.01, 0.001)), list(theta = c(0.001, 0.001)),
          list(theta = c(0.005, 0.005)), list(theta = c(0.01, 0.01)))
      } else {
        stan_inits <- list(
          list(theta = 0.01), list(theta = 0.1),
          list(theta = 0.001), list(theta = 0.005))
        
      }
    }
    if (k == 0) {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = if(!missing(slope)) halfnormal(slope[1], slope[2]),
        a1 = a1, P1 = as.matrix(P1), distribution = "poisson")
      ## refine initial mode
      model$initial_mode[] <- smoother(model)$alphahat[, 1]
      stan_data <- list(y = y, n = length(y),
        a1 = structure(a1, dim = length(a1)), P1 = as.matrix(P1), 
        sd_prior_means = 
          c(level[1], if(!missing(slope)) slope[1]), dim = 1 + !missing(slope),
        sd_prior_sds = 
          c(level[2], if(!missing(slope)) slope[2]), dim = 1 + !missing(slope),
        initial_mode = c(model$initial_mode, -1e300), distribution = 1L, 
        max_iter = 100, conv_tol = 1e-8)
      
      fit <- sampling(
        if(missing(slope)) stanmodels$ll_approx else stanmodels$llt_approx,
        data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, cores = n_threads, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = if(!missing(slope)) halfnormal(slope[1],slope[2]),
        a1 = a1, P1 = as.matrix(P1), distribution = "poisson", xreg = xreg,
        beta = normal(0, beta[, 1], beta[, 2]))
      ## refine initial mode
      model$initial_mode[] <- smoother(model)$alphahat[, 1]
      
      stan_data <- list(y = y, n = length(y), k = k,
        a1 = structure(a1, dim = length(a1)), P1 = as.matrix(P1), 
        sd_prior_means = 
          c(level[1], if(!missing(slope)) slope[1]), dim = 1 + !missing(slope),
        sd_prior_sds = c(level[2], if(!missing(slope)) slope[2]), dim = 1 + !missing(slope),
        beta_prior_means = beta[, 1], beta_prior_sds = beta[, 2],
        initial_mode = c(model$initial_mode, -1e300), xreg = xreg, distribution = 1L,
        max_iter = 100, conv_tol = 1e-8)
      if (k == 1) {
        dim(stan_data$beta_prior_means) <- dim(stan_data$beta_prior_sds) <- 1
      }
      fit <- sampling( 
        if(missing(slope)) stanmodels$x_ll_approx else stanmodels$x_llt_approx,
        data = stan_data,
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
        sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = a1, P1 = as.matrix(P1),
        distribution = "poisson")
      ## refine initial mode
      model$initial_mode[] <- smoother(model)$alphahat[, 1] + smoother(model)$alphahat[, 3]
      
      stan_data <- list(y = y, n = length(y), period = period,
        a1 = structure(a1, dim = length(a1)), P1 = as.matrix(P1), sd_prior_means = c(level[1], slope[1], seasonal[1]),
        sd_prior_sds = c(level[2], slope[2], seasonal[2]),
        initial_mode = c(model$initial_mode, -1e300), distribution = 1L,
        max_iter = 100, conv_tol = 1e-8)
      
      fit <- sampling(stanmodels$bsm_approx, data = stan_data,
        chains = length(stan_inits), init = stan_inits,
        pars = "Rt", include = FALSE,
        iter = iter,  thin = thin, cores = n_threads, ...)
      
    } else {
      model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
        sd_slope = halfnormal(slope[1], slope[2]),
        sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = a1, P1 = as.matrix(P1),
        distribution = "poisson", xreg = xreg,
        beta = normal(0, beta[, 1], beta[, 2]))
      
      ## refine initial mode
      model$initial_mode[] <- smoother(model)$alphahat[, 1] + smoother(model)$alphahat[, 3]
      
      stan_data <- list(y = y, n = length(y), k = k, period = period,
        a1 = structure(a1, dim = length(a1)), P1 = as.matrix(P1), 
        sd_prior_means = c(level[1], slope[1], seasonal[1]),
        sd_prior_sds = c(level[2], slope[2], seasonal[2]),
        initial_mode = c(model$initial_mode, -1e300), distribution = 1L,
        xreg = xreg, beta_prior_means = beta[, 1], beta_prior_sds = beta[, 2], 
        max_iter = 100, conv_tol = 1e-8)
      if (k == 1) {
        dim(stan_data$beta_prior_means) <- dim(stan_data$beta_prior_sds) <- 1
      }
      fit <- sampling(stanmodels$x_bsm_approx, data = stan_data,
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
      select_(stan_out, "lp__")[, 1], select_(stan_out, "jacobian")[, 1], 
      nsim_states, n_threads, is_seed)
    c_time <- proc.time() - c_time
  } else {
    
    stan_out <- as.data.frame(apply( 
      extract(fit, c("theta", "beta", "lp__", "jacobian", "approx_results"), permuted = FALSE), 
      3, identity))
    
    c_time <- proc.time()
    correction <- is_correction(model, as.matrix(select(stan_out, starts_with("approx_results"))),
      as.matrix(cbind(select(stan_out, starts_with("theta")), select(stan_out, starts_with("beta")))), 
      select_(stan_out, "lp__")[, 1], select_(stan_out, "jacobian")[, 1], 
      nsim_states, n_threads, is_seed)
    c_time <- proc.time() - c_time
  }
  
  w <- exp(c(correction$weights))
  n_chains <- length(stan_inits)
  N <- length(w) / n_chains
  
  
  theta <- as.matrix(select(stan_out, starts_with("theta")))
  mean_theta <- weighted_mean(theta, w)
  ess_ar <- matrix(NA, length(mean_theta), n_chains)
  for(i in 1:n_chains) {
    ind <- (i-1)*N + 1:N
    mean_theta_i <-  weighted_mean(theta[ind,], w[ind])
    sd_theta_i <- if(length(mean_theta) > 1) sqrt(diag(weighted_var(theta[ind,], w[ind]))) else sqrt(weighted_var(theta[ind,], w[ind]))
    spec_i <- sapply(1:length(sd_theta_i), function(j) 
      spectrum0.ar((theta[ind, j] - mean_theta_i[j]) * w[ind])$spec)
    se_theta_i <- sqrt(spec_i / length(w[ind])) / mean(w[ind])
    ess_ar[, i] <- (sd_theta_i / se_theta_i)^2
    
  }
  ess_theta <- rowSums(ess_ar)
  if(k > 0) {
    beta <- as.matrix(select(stan_out, starts_with("beta")))
    mean_beta <- weighted_mean(beta, w)
    ess_ar <- matrix(NA, length(mean_beta), n_chains)
    for(i in 1:n_chains) {
      ind <- (i-1)*N + 1:N
      mean_beta_i <-  weighted_mean(beta[ind,], w[ind])
      sd_beta_i <- if(k > 1) sqrt(diag(weighted_var(beta[ind,], w[ind]))) else sqrt(weighted_var(beta[ind,], w[ind]))
      spec_i <- sapply(1:length(sd_beta_i), function(j) 
        spectrum0.ar((beta[ind, j] - mean_beta_i[j]) * w[ind])$spec)
      se_beta_i <- sqrt(spec_i / length(w[ind])) / mean(w[ind])
      ess_ar[, i] <- (sd_beta_i / se_beta_i)^2
      
    }
    ess_beta <- rowSums(ess_ar)
  } else {
    beta <- NULL
    ess_beta <- NULL
    mean_beta <- NULL
  }
  
  if (simplify_output) {
    list(theta = theta, mean_theta = mean_theta, ess_theta = ess_theta, 
      beta = beta, mean_beta = mean_beta, ess_beta = ess_beta,
      states = correction$alpha, weights = w,
      stan_time = sum(get_elapsed_time(fit)), correction_time = c_time[3])
  } else {
    list(stan_fit = fit, theta = as.matrix(select(stan_out, starts_with("theta"))), 
      beta = if(k > 0) as.matrix(select(stan_out, starts_with("beta"))),
      correction = correction, stan_time = get_elapsed_time(fit), correction_time = c_time)
  }
}

