#include "base.h"
#include "is_correction_bsf.h"
#include "ung_bsm.h"

// [[Rcpp::export]]
Rcpp::List R_is_correction_bsf(const Rcpp::List& model_,
  const arma::mat& theta, const arma::vec& approx_posterior,
  const unsigned int nsim_states, const arma::vec& prior, 
  const unsigned int n_threads, const unsigned int seed) {
  
  ung_bsm model(clone(model_), seed);
  return is_correction_bsf(model, theta, approx_posterior, nsim_states, 
    prior, n_threads);
}
