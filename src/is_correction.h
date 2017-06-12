#ifndef IS_H
#define IS_H

#include <RcppArmadillo.h>
class ung_bsm;

Rcpp::List is_correction_global(ung_bsm model, 
  const arma::mat& approx_y, const arma::mat& approx_var_y, 
  const arma::mat& theta, 
  const arma::vec& approx_posterior, const arma::vec& jacobian, 
  const unsigned int nsim_states, const unsigned int n_threads);

Rcpp::List is_correction(ung_bsm model, 
  const arma::mat& approx_y, const arma::mat& approx_var_y, 
  const arma::mat& scales, const arma::mat& theta, 
  const arma::vec& approx_posterior, const arma::vec& jacobian, 
  const unsigned int nsim_states, const unsigned int n_threads);

void state_sampler(ung_bsm model, const unsigned int nsim_states, 
  const arma::mat& theta, arma::cube& alpha, arma::vec& weights, 
  const arma::mat& approx_y, const arma::mat& approx_var_y, const arma::mat& scales);


void state_sampler_global(ung_bsm model, const unsigned int nsim_states, 
  const arma::mat& theta, arma::cube& alpha, arma::vec& weights, 
  const arma::vec& approx_y, const arma::vec& approx_var_y);
#endif
