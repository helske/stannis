#' @export
is_correction <- function(object, theta, posterior, prior, nsim_states, n_threads = 1, seed = 1) {
  object$distribution <- pmatch(object$distribution, c("poisson", "binomial", "negative binomial"))
  R_is_correction_bsf(object, t(theta), posterior, nsim_states, prior, n_threads, seed)
}
