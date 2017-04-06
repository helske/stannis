// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// conditional_cov
void conditional_cov(arma::cube& Vt, arma::cube& Ct, const bool use_svd);
RcppExport SEXP stannis_conditional_cov(SEXP VtSEXP, SEXP CtSEXP, SEXP use_svdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type Vt(VtSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Ct(CtSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_svd(use_svdSEXP);
    conditional_cov(Vt, Ct, use_svd);
    return R_NilValue;
END_RCPP
}
// R_is_correction_bsf
Rcpp::List R_is_correction_bsf(const Rcpp::List& model_, const arma::mat& theta, const arma::vec& approx_posterior, const unsigned int nsim_states, const arma::vec& prior, const unsigned int n_threads, const unsigned int seed);
RcppExport SEXP stannis_R_is_correction_bsf(SEXP model_SEXP, SEXP thetaSEXP, SEXP approx_posteriorSEXP, SEXP nsim_statesSEXP, SEXP priorSEXP, SEXP n_threadsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type approx_posterior(approx_posteriorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type nsim_states(nsim_statesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(R_is_correction_bsf(model_, theta, approx_posterior, nsim_states, prior, n_threads, seed));
    return rcpp_result_gen;
END_RCPP
}
// R_is_correction_psi
Rcpp::List R_is_correction_psi(const Rcpp::List& model_, const arma::mat& approx_y, const arma::mat approx_var_y, const arma::mat scales, const arma::mat& theta, const arma::vec& approx_posterior, const arma::vec& jacobian, const unsigned int nsim_states, const unsigned int n_threads, const unsigned int seed);
RcppExport SEXP stannis_R_is_correction_psi(SEXP model_SEXP, SEXP approx_ySEXP, SEXP approx_var_ySEXP, SEXP scalesSEXP, SEXP thetaSEXP, SEXP approx_posteriorSEXP, SEXP jacobianSEXP, SEXP nsim_statesSEXP, SEXP n_threadsSEXP, SEXP seedSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(R_is_correction_psi(model_, approx_y, approx_var_y, scales, theta, approx_posterior, jacobian, nsim_states, n_threads, seed));
    return rcpp_result_gen;
END_RCPP
}
