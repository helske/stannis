#include "ung_bsm.h"
#include "is_correction.h"
#include "filter_smoother.h"

Rcpp::List is_correction(ung_bsm model, 
  const arma::mat& approx_y, const arma::mat& approx_var_y, const arma::mat& scales, 
  const arma::mat& theta, const arma::vec& approx_posterior, 
  const arma::vec& jacobian, const unsigned int nsim_states, 
  const unsigned int n_threads) {
  
  unsigned int n_samples = approx_posterior.n_elem;
  unsigned int m = model.m;
  unsigned int n = model.n;
  arma::cube alpha(n, m, n_samples);
  arma::vec weights(n_samples);
  
  if(n_threads > 1) {
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) default(none) \
    shared(alpha, weights, n_samples, theta, approx_y, approx_var_y, scales) firstprivate(model)
    {
      model.engine = std::mt19937(omp_get_thread_num() + 1);
      unsigned thread_size = floor(n_samples / n_threads);
      unsigned int start = omp_get_thread_num() * thread_size;
      unsigned int end = (omp_get_thread_num() + 1) * thread_size - 1;
      if(omp_get_thread_num() == (n_threads - 1)) {
        end = n_samples - 1;
      }
      
      arma::mat theta_piece = theta(arma::span::all, arma::span(start, end));
      arma::cube alpha_piece(model.n, model.m, thread_size);
      arma::vec weights_piece(thread_size);
      arma::mat approx_y_piece = approx_y(arma::span::all, arma::span(start, end));
      arma::mat approx_var_y_piece = approx_var_y(arma::span::all, arma::span(start, end));
      arma::mat scales_piece = scales(arma::span::all, arma::span(start,end));
      state_sampler(model, nsim_states, theta_piece, 
        alpha_piece, weights_piece, approx_y_piece, approx_var_y_piece, scales_piece);
      
      alpha.slices(start, end) = alpha_piece;
      weights.subvec(start, end) = weights_piece;
    }
#else
    state_sampler(model, nsim_states, theta, alpha, 
      weights, approx_y, approx_var_y, scales);
#endif
  } else {
    state_sampler(model, nsim_states, theta, alpha, 
      weights, approx_y, approx_var_y, scales);
  }
  
  arma::vec posterior = approx_posterior + weights - jacobian;
  //weights -= jacobian;
  // double wmax = weights.max();
  // weights = exp(weights - wmax);
  // double wsum = arma::sum(weights);
  // weights /= wsum; 
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
    Rcpp::Named("weights") = weights, Rcpp::Named("posterior") = posterior);
}

void state_sampler(ung_bsm model, const unsigned int nsim_states, 
  const arma::mat& theta, arma::cube& alpha, arma::vec& weights, 
  const arma::mat& approx_y, const arma::mat& approx_var_y, const arma::mat& scales) {
  
  for (unsigned int i = 0; i < theta.n_cols; i++) {
    
    arma::vec theta_i = theta.col(i);
    model.set_theta(theta_i);
    arma::cube alpha_i(model.m, model.n, nsim_states);
    arma::mat weights_i(nsim_states, model.n);
    arma::umat indices(nsim_states, model.n - 1);
    weights(i) = model.psi_filter(approx_y.col(i), approx_var_y.col(i), 
      scales.col(i), nsim_states, alpha_i, weights_i, indices);
    filter_smoother(alpha_i, indices);
    arma::vec w = weights_i.col(model.n - 1);
    std::discrete_distribution<> sample(w.begin(), w.end());
    alpha.slice(i) = alpha_i.slice(sample(model.engine)).t();
  }
}



Rcpp::List is_correction_global(ung_bsm model, 
  const arma::mat& approx_y, const arma::mat& approx_var_y,
  const arma::mat& theta, const arma::vec& approx_posterior, 
  const arma::vec& jacobian, const unsigned int nsim_states, 
  const unsigned int n_threads) {
  
  unsigned int n_samples = approx_posterior.n_elem;
  unsigned int m = model.m;
  unsigned int n = model.n;
  unsigned int n_chains = approx_y.n_cols;
  arma::cube alpha(n, m, n_samples);
  arma::vec weights(n_samples);
  
  if(n_chains > 1) {
#ifdef _OPENMP
#pragma omp parallel num_threads(n_chains) default(none) \
    shared(alpha, weights, n_samples, theta, approx_y, approx_var_y, n_chains) firstprivate(model)
    {
      model.engine = std::mt19937(omp_get_thread_num() + 1);
      unsigned thread_size = floor(n_samples / n_chains);
      unsigned int start = omp_get_thread_num() * thread_size;
      unsigned int end = (omp_get_thread_num() + 1) * thread_size - 1;
      if(omp_get_thread_num() == (n_chains - 1)) {
        end = n_samples - 1;
      }
      
      arma::mat theta_piece = theta(arma::span::all, arma::span(start, end));
      arma::cube alpha_piece(model.n, model.m, thread_size);
      arma::vec weights_piece(thread_size);
      
      state_sampler_global(model, nsim_states, theta_piece, 
        alpha_piece, weights_piece, approx_y.col(omp_get_thread_num()), approx_var_y.col(omp_get_thread_num()));
      
      alpha.slices(start, end) = alpha_piece;
      weights.subvec(start, end) = weights_piece;
    }
#else
    for(unsigned int i = 0; i < n_chains; i++) {
      unsigned thread_size = std::floor(n_samples / n_chains);
      unsigned int start = i * thread_size;
      unsigned int end = (i + 1) * thread_size - 1;
     
      arma::mat theta_piece = theta(arma::span::all, arma::span(start, end));
      arma::cube alpha_piece(model.n, model.m, thread_size);
      arma::vec weights_piece(thread_size);
      
      state_sampler_global(model, nsim_states, theta_piece, 
        alpha_piece, weights_piece, approx_y.col(i), approx_var_y.col(i));
      
      alpha.slices(start, end) = alpha_piece;
      weights.subvec(start, end) = weights_piece;
    }
#endif
  } else {
    state_sampler_global(model, nsim_states, theta, alpha, 
      weights, approx_y, approx_var_y);
  }
  
  arma::vec posterior = approx_posterior + weights - jacobian;
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
    Rcpp::Named("weights") = weights, Rcpp::Named("posterior") = posterior);
}

void state_sampler_global(ung_bsm model, const unsigned int nsim_states, 
  const arma::mat& theta, arma::cube& alpha, arma::vec& weights, 
  const arma::vec& approx_y, const arma::vec& approx_var_y) {
  
  for (unsigned int i = 0; i < theta.n_cols; i++) {
    
    arma::vec theta_i = theta.col(i);
    model.set_theta(theta_i);
    arma::cube alpha_i(model.m, model.n, nsim_states);
    arma::mat weights_i(nsim_states, model.n);
    arma::umat indices(nsim_states, model.n - 1);
    weights(i) = model.psi_filter(approx_y, approx_var_y, 
      arma::vec(model.n, arma::fill::zeros), nsim_states, alpha_i, weights_i, indices);
    filter_smoother(alpha_i, indices);
    arma::vec w = weights_i.col(model.n - 1);
    std::discrete_distribution<> sample(w.begin(), w.end());
    alpha.slice(i) = alpha_i.slice(sample(model.engine)).t();
  }
}
