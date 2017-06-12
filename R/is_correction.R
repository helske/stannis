#' @export
is_correction <- function(object, approx_model, theta, posterior, jacobian,
  nsim_states, n_threads = 1, seed = 1) {
  object$distribution <- pmatch(object$distribution, c("poisson", "binomial", "negative binomial"))
  n <- length(object$y)
  R_is_correction(object, t(approx_model[,1:n]), t(approx_model[,n + 1:n]),
    t(approx_model[, 2 * n + 1:n]), t(theta), posterior, jacobian, nsim_states, n_threads, seed)
}

#' @export
is_correction_global <- function(object, approx_y, approx_y_var, theta, posterior, jacobian,
  nsim_states, n_threads = 1, seed = 1) {
  object$distribution <- pmatch(object$distribution, c("poisson", "binomial", "negative binomial"))
  n <- length(object$y)
  R_is_correction_global(object, approx_y, approx_y_var, t(theta), posterior, jacobian, 
    nsim_states, n_threads, seed)
}

