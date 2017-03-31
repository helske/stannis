#ifndef IS_BSF_H
#define IS_BSF_H

#include <RcppArmadillo.h>
class ung_bsm;

Rcpp::List is_correction_bsf(ung_bsm model, const arma::mat& theta, 
  const arma::vec& approx_posterior, const unsigned int nsim_states, 
  const arma::vec& prior, const unsigned int n_threads);

void state_sampler_bsf(ung_bsm model, const unsigned int nsim_states, 
  const arma::mat& theta, const arma::vec& approx_posterior, 
  const arma::vec& prior, arma::cube& alpha, arma::vec& weights);

#endif
