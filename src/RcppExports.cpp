// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// R_is_correction
Rcpp::List R_is_correction(const Rcpp::List& model_, const arma::mat& approx_y, const arma::mat approx_var_y, const arma::mat scales, const arma::mat& theta, const arma::vec& approx_posterior, const arma::vec& jacobian, const unsigned int nsim_states, const unsigned int n_threads, const unsigned int seed);
RcppExport SEXP stannis_R_is_correction(SEXP model_SEXP, SEXP approx_ySEXP, SEXP approx_var_ySEXP, SEXP scalesSEXP, SEXP thetaSEXP, SEXP approx_posteriorSEXP, SEXP jacobianSEXP, SEXP nsim_statesSEXP, SEXP n_threadsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type approx_y(approx_ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type approx_var_y(approx_var_ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type approx_posterior(approx_posteriorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type jacobian(jacobianSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type nsim_states(nsim_statesSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(R_is_correction(model_, approx_y, approx_var_y, scales, theta, approx_posterior, jacobian, nsim_states, n_threads, seed));
    return rcpp_result_gen;
END_RCPP
}
