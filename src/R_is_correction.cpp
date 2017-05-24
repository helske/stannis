#include "base.h"
#include "is_correction.h"
#include "ung_bsm.h"

//[[Rcpp::export]]
Rcpp::List is_correction(const Rcpp::List& model_,
  const arma::mat& approx_y, const arma::mat approx_var_y, const arma::mat scales,
  const arma::mat& theta, const arma::vec& approx_posterior,
  const arma::vec& jacobian, const unsigned int nsim_states,
  const unsigned int n_threads, const unsigned int seed) {

  ung_bsm model(clone(model_), seed);
  return is_correction(model, approx_y, approx_var_y, scales, theta, approx_posterior,
    jacobian, nsim_states, n_threads);
}
