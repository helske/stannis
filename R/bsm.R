#' #' @export
#' bsm <- function(formula, data, level, slope, seasonal,
#'                 distribution, x1, P1, beta,
#'                 iter = 7500, thin = 5,
#'                 nsim_states = 10, stan_inits, ...) {
#' 
#'   distr <- pmatch(distribution, c("poisson", "binomial", "negative binomial"))
#' 
#'   # build xreg
#'   if(inherits(formula, "formula")) {
#'     if (missing(data)) {
#'       data <- environment(formula)
#'       period <- NULL
#'     } else {
#'       period <- tsp(data)[3]
#'       data <- as.data.frame(data)
#'     }
#'     mf <- match.call(expand.dots = FALSE)
#'     mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
#'     mf$drop.unused.levels <- TRUE
#'     mf[[1L]] <- quote(stats::model.frame)
#'     mf$na.action <- as.name("na.pass")
#'     mf <- eval(mf, parent.frame())
#'     # remove the intercept
#'     attr(attr(mf, "terms"), "intercept") <- 0L
#'     y <- model.response(mf, "numeric")
#'     xreg <- model.matrix(attr(mf, "terms"), mf)
#'     k <- ncol(xreg)
#'     if (is.null(period)) {
#'       period <- frequency(y)
#'     }
#'   } else {
#'     y <- formula
#'     k <- 0L
#'     period <- frequency(y)
#'   }
#' 
#'   ## local level or local linear trend model
#' 
#'   if (period == 1) {
#' 
#'     ## initial values
#'     if (missing(stan_inits)) {
#'       if (!missing(slope)) {
#'         stan_inits <- list(
#'           list(theta = c(0.01, 0.0001)), list(theta = c(0.1, 0.1)),
#'           list(theta = c(0.0001, 0.01)), list(theta = c(0.001, 0.001)))
#'       } else {
#'         stan_inits <- list(
#'           list(theta = c(0.01)), list(theta = c(0.1)),
#'           list(theta = c(0.0001)), list(theta = c(0.001)))
#'       }
#'     }
#'     if (k == 0) {
#'       model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
#'                       sd_slope = if(!missing(slope)) halfnormal(slope[1],slope[2]),
#'                       a1 = x1, P1 = P1, distribution = distribution)
#' 
#'       stan_data <- list(y = y, n = length(y),
#'                         x1 = x1, P1 = P1, sd_prior_means = c(level[1], if(!missing(slope)) slope[1]),
#'                         sd_prior_sds = c(level[2], if(!missing(slope)) slope[2]),
#'                         initial_mode = c(model$initial_mode, -1e300), distribution = distr)
#' 
#'       fit <- sampling(
#'         if(missing(slope)) stannis:::stanmodels$ll_approx else stannis:::stanmodels$llt_approx,
#'            data = stan_data,
#'            chains = length(stan_inits), init = stan_inits,
#'            pars = "Rt", include = FALSE,
#'            iter = iter,  thin = thin, ...)
#' 
#'     } else {
#'       model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
#'                       sd_slope = if(!missing(slope)) halfnormal(slope[1],slope[2]),
#'                       a1 = x1, P1 = P1, distribution = distribution, xreg = xreg,
#'                       beta = normal(0, beta[,1], beta[,2]))
#' 
#'       stan_data <- list(y = y, n = length(y), k = k,
#'                         x1 = x1, P1 = P1, sd_prior_means = c(level[1], if(!missing(slope)) slope[1]),
#'                         sd_prior_sds = c(level[2], if(!missing(slope)) slope[2]),
#'                         beta_prior_means = beta[,1], beta_prior_sds = beta[,2],
#'                         initial_mode = c(model$initial_mode, -1e300), xreg = xreg, distribution = distr)
#' 
#'       fit <- sampling(stannis:::stanmodels$x_llt_approx, data = stan_data,
#'                       chains = length(stan_inits), init = stan_inits,
#'                       pars = c("Rt", "xbeta"), include = FALSE,
#'                       iter = iter,  thin = thin, ...)
#'     }
#' 
#'   } else {
#'     # BSM, note we always include slope
#'     if (missing(stan_inits)) {
#'         stan_inits <- list(
#'           list(theta = rep(0.001, 3)), list(theta = c(0.01, 0.001, 0.001)),
#'           list(theta = c(0.001, 0.001, 0.1)),list(theta = c(0.01, 0.01, 0.001)))
#'     }
#'     if(k == 0L) {
#'       model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
#'                       sd_slope = halfnormal(slope[1], slope[2]),
#'                       sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = x1, P1 = P1,
#'                       distribution = distribution)
#' 
#'       stan_data <- list(y = y, n = length(y), period = period,
#'                         x1 = x1, P1 = P1, sd_prior_means = c(level[1], slope[1], seasonal[1]),
#'                         sd_prior_sds = c(level[2], slope[2], seasonal[2]),
#'                         initial_mode = c(model$initial_mode, -1e300), distribution = distr)
#' 
#'       fit <- sampling(stannis:::stanmodels$bsm_approx, data = stan_data,
#'                       chains = length(stan_inits), init = stan_inits,
#'                       pars = "Rt", include = FALSE,
#'                       iter = iter,  thin = thin, ...)
#' 
#'     } else {
#'       model <- ng_bsm(y, sd_level = halfnormal(level[1], level[2]),
#'                       sd_slope = halfnormal(slope[1], slope[2]),
#'                       sd_seasonal = halfnormal(seasonal[1], seasonal[2]),  a1 = x1, P1 = P1,
#'                       distribution = distribution, xreg = xreg,
#'                       beta = normal(0, beta[, 1], beta[, 2]))
#' 
#' 
#' 
#'       stan_data <- list(y = y, n = length(y), k = k, period = period,
#'                         x1 = x1, P1 = P1, sd_prior_means = c(level[1], slope[1], seasonal[1]),
#'                         sd_prior_sds = c(level[2], slope[2], seasonal[2]),
#'                         initial_mode = c(model$initial_mode, -1e300), distribution = distr,
#'                         xreg = xreg, beta_prior_means = beta[, 1], beta_prior_sds = beta[, 2])
#' 
#'       fit <- sampling(stannis:::stanmodels$x_bsm_approx, data = stan_data,
#'                       chains = length(stan_inits), init = stan_inits,
#'                       pars = c("Rt", "xbeta"), include = FALSE,
#'                       iter = iter, thin = thin, ...)
#'     }
#'   }
#' 
#' 
#'   if (k==0) {
#'     c_time <- proc.time()
#'     stan_out <- extract(fit, c("lp__", "theta", "jacobian", "approx_results"))
#'     correction <- is_correction(model, stan_out$approx_results,
#'                                     stan_out$theta, stan_out$lp__, stan_out$jacobian, nsim_states)
#'     c_time <- proc.time() - c_time
#'   } else {
#'     c_time <- proc.time()
#'     stan_out <- extract(fit, c("lp__", "theta", "beta", "jacobian", "approx_results"))
#'     correction <- is_correction(model, stan_out$approx_results,
#'                                     cbind(stan_out$theta, stan_out$beta), stan_out$lp__, stan_out$jacobian, nsim_states)
#'     c_time <- proc.time() - c_time
#'   }
#'   list(stan_fit = fit, correction = correction, s_time = get_elapsed_time(fit), c_time = c_time)
#' }
#' 
