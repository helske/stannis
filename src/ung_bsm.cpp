#include "ung_bsm.h"
#include "stratified_sample.h"
#include "distr_consts.h"
#include "conditional_dist.h"

ung_bsm::ung_bsm(const Rcpp::List& model, const unsigned int seed) :
  y(Rcpp::as<arma::vec>(model["y"])), Z(Rcpp::as<arma::mat>(model["Z"])),
  T(Rcpp::as<arma::cube>(model["T"])), R(Rcpp::as<arma::cube>(model["R"])), 
  a1(Rcpp::as<arma::vec>(model["a1"])), P1(Rcpp::as<arma::mat>(model["P1"])), 
  xreg(Rcpp::as<arma::mat>(model["xreg"])), beta(Rcpp::as<arma::vec>(model["coefs"])), 
  D(Rcpp::as<arma::vec>(model["obs_intercept"])), C(Rcpp::as<arma::mat>(model["state_intercept"])),
  Ztv(Z.n_cols > 1), Ttv(T.n_slices > 1), Rtv(R.n_slices > 1), Dtv(D.n_elem > 1), 
  Ctv(C.n_cols > 1),
  n(y.n_elem), m(a1.n_elem), k(R.n_cols), RR(arma::cube(m, m, Rtv * (n - 1) + 1)), 
  xbeta(arma::vec(n, arma::fill::zeros)), engine(seed),
  phi(model["phi"]), 
  u(Rcpp::as<arma::vec>(model["u"])), distribution(model["distribution"]), 
  phi_est(Rcpp::as<bool>(model["phi_est"])),
  seed(seed), slope(Rcpp::as<bool>(model["slope"])),
  seasonal(Rcpp::as<bool>(model["seasonal"])),
  noise(Rcpp::as<bool>(model["noise"])),
  fixed(Rcpp::as<arma::uvec>(model["fixed"])), level_est(fixed(0) == 0),
  slope_est(slope && fixed(1) == 0), seasonal_est(seasonal && fixed(2) == 0) {
  
  if(xreg.n_cols > 0) {
    compute_xbeta();
  }
  compute_RR();
}

void ung_bsm::compute_RR(){
  for (unsigned int t = 0; t < R.n_slices; t++) {
    RR.slice(t) = R.slice(t * Rtv) * R.slice(t * Rtv).t();
  }
}

void ung_bsm::set_theta(const arma::vec& theta) {
  
  if (sum(fixed) < 3 || noise || phi_est) {
    
    // sd_level
    if (level_est) {
      R(0, 0, 0) = theta(0);
    }
    // sd_slope
    if (slope_est) {
      R(1, 1, 0) = theta(level_est);
    }
    // sd_seasonal
    if (seasonal_est) {
      R(1 + slope, 1 + slope, 0) =
        theta(level_est + slope_est);
    }
    if(noise) {
      R(m - 1, 1 + slope + seasonal, 0) =
        theta(level_est + slope_est + seasonal_est);
      P1(m - 1, m - 1) = std::pow(theta(level_est + slope_est + seasonal_est), 2);
    }
    compute_RR();
  }
  if(phi_est) {
    phi = theta(level_est + slope_est + seasonal_est + noise);
  }
  
  if(xreg.n_cols > 0) {
    beta = theta.subvec(theta.n_elem - xreg.n_cols, theta.n_elem - 1);
    compute_xbeta();
  }
}

arma::vec ung_bsm::get_theta(void) const {
  
  unsigned int npar = level_est + slope_est + seasonal_est + noise +
    xreg.n_cols + phi_est;
  
  arma::vec theta(npar);
  
  if (sum(fixed) < 3 || noise) {
    // sd_level
    if (level_est) {
      theta(0) = R(0, 0, 0);
    }
    // sd_slope
    if (slope_est) {
      theta(level_est) = R(1, 1, 0);
    }
    // sd_seasonal
    if (seasonal_est) {
      theta(level_est + slope_est) =
        R(1 + slope, 1 + slope, 0);
    }
    if (noise) {
      theta(level_est + slope_est + seasonal_est) =
        R(m - 1, 1 + slope + seasonal, 0);
    }
  }
  
  if(phi_est) {
    theta(level_est + slope_est + seasonal_est + noise) = phi;
  }
  
  if(xreg.n_cols > 0) {
    theta.subvec(theta.n_elem - xreg.n_cols, theta.n_elem - 1) = beta;
  }
  return theta;
}


// Logarithms of _unnormalized_ densities g(y_t | alpha_t)
/*
 * t:             Time point where the densities are computed
 * alpha:         Simulated particles
 */
arma::vec ung_bsm::log_obs_density(const unsigned int t, 
  const arma::cube& alpha) const {
  
  arma::vec weights(alpha.n_slices, arma::fill::zeros);
  
  switch(distribution) {
  case 0  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double simsignal = alpha(0, t, i);
      weights(i) = -0.5 * (simsignal + pow(y(t) / phi, 2) * exp(-simsignal));
    }
    break;
  case 1  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double simsignal = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t) + xbeta(t));
      weights(i) = y(t) * simsignal  - u(t) * exp(simsignal);
    }
    break;
  case 2  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double simsignal = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t) + xbeta(t));
      weights(i) = y(t) * simsignal - u(t) * log1p(exp(simsignal));
    }
    break;
  case 3  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double simsignal = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t) + xbeta(t));
      weights(i) = y(t) * simsignal - (y(t) + phi) * 
        log(phi + u(t) * exp(simsignal));
    }
    break;
    
  }
  return weights;
}

double ung_bsm::bsf_filter(const unsigned int nsim, arma::cube& alpha,
  arma::mat& weights, arma::umat& indices) {
  
  // arma::mat U(m, m);
  // arma::mat V(m, m);
  // arma::vec s(m);
  // arma::svd_econ(U, s, V, P1, "left");
  // arma::uvec nonzero = arma::find(s > (arma::datum::eps * m * s(0)));
  // arma::mat L = arma::diagmat(1.0 / s(nonzero)) U
  arma::uvec nonzero = arma::find(P1.diag() > 0);
  arma::mat L_P1(m, m, arma::fill::zeros);
  if (nonzero.n_elem > 0) {
    L_P1.submat(nonzero, nonzero) =
      arma::chol(P1.submat(nonzero, nonzero), "lower");
  }
  std::normal_distribution<> normal(0.0, 1.0);
  for (unsigned int i = 0; i < nsim; i++) {
    arma::vec um(m);
    for(unsigned int j = 0; j < m; j++) {
      um(j) = normal(engine);
    }
    alpha.slice(i).col(0) = a1 + L_P1 * um;
  }
  
  std::uniform_real_distribution<> unif(0.0, 1.0);
  arma::vec normalized_weights(nsim);
  double loglik = 0.0;
  
  
  weights.col(0) = log_obs_density(0, alpha);
  double max_weight = weights.col(0).max();
  weights.col(0) = exp(weights.col(0) - max_weight);
  double sum_weights = arma::sum(weights.col(0));
  if(sum_weights > 0.0){
    normalized_weights = weights.col(0) / sum_weights;
  } else {
    return -arma::datum::inf;
  }
  loglik = max_weight + log(sum_weights / nsim);
  
  for (unsigned int t = 0; t < (n - 1); t++) {
    
    arma::vec r(nsim);
    for (unsigned int i = 0; i < nsim; i++) {
      r(i) = unif(engine);
    }
    
    indices.col(t) = stratified_sample(normalized_weights, r, nsim);
    
    arma::mat alphatmp(m, nsim);
    
    for (unsigned int i = 0; i < nsim; i++) {
      alphatmp.col(i) = alpha.slice(indices(i, t)).col(t);
    }
    
    for (unsigned int i = 0; i < nsim; i++) {
      arma::vec uk(k);
      for(unsigned int j = 0; j < k; j++) {
        uk(j) = normal(engine);
      }
      alpha.slice(i).col(t + 1) = C.col(t * Ctv) + 
        T.slice(t * Ttv) * alphatmp.col(i) + R.slice(t * Rtv) * uk;
    }
    
    weights.col(t + 1) = log_obs_density(t + 1, alpha);
    
    double max_weight = weights.col(t + 1).max();
    weights.col(t + 1) = exp(weights.col(t + 1) - max_weight);
    double sum_weights = arma::sum(weights.col(t + 1));
    if(sum_weights > 0.0){
      normalized_weights = weights.col(t + 1) / sum_weights;
    } else {
      return -arma::datum::inf;
    }
    loglik += max_weight + log(sum_weights / nsim);
    
  }
  // constant part of the log-likelihood
  switch(distribution) {
  case 0 :
    loglik += y.n_elem * norm_log_const(phi);
    break;
  case 1 : {
      loglik += poisson_log_const(y, u);
    } break;
  case 2 : {
    loglik += binomial_log_const(y, u);
  } break;
  case 3 : {
    loglik += negbin_log_const(y, u, phi);
  } break;
  }
  return loglik;
}

arma::vec ung_bsm::log_weights(const arma::vec& approx_y, const arma::vec& HH, 
  const unsigned int t, const arma::cube& alpha) const {
  
  arma::vec weights(alpha.n_slices, arma::fill::zeros);
  
  switch(distribution) {
  case 0  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double zt = alpha(0, t, i);
      weights(i) = -0.5 * (zt + pow(y(t) / phi, 2) * exp(-zt)) +
        0.5 * std::pow(approx_y(t) - zt, 2) / HH(t) + log(HH(t));
    }
    break;
  case 1  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double zt = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t));
      weights(i) = y(t) * (zt + xbeta(t))  - u(t) * exp(zt + xbeta(t)) +
        0.5 * std::pow(approx_y(t) - zt, 2) / HH(t) + log(HH(t));
    }
    break;
  case 2  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double zt = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t));
      weights(i) = y(t) * (zt + xbeta(t)) - u(t) * log1p(exp(zt + xbeta(t))) +
        0.5 * std::pow(approx_y(t) - zt, 2) / HH(t) + log(HH(t));
    }
    break;
  case 3  :
    for (unsigned int i = 0; i < alpha.n_slices; i++) {
      double zt = arma::as_scalar(Z.col(t * Ztv).t() *
        alpha.slice(i).col(t));
      weights(i) = y(t) * (zt + xbeta(t)) - (y(t) + phi) *
        log(phi + u(t) * exp(zt + xbeta(t))) +
        0.5 * std::pow(approx_y(t) - zt, 2) / HH(t) + log(HH(t));
    }
    break;
  }
  
  return weights;
}


// approximate smoother which returns also cov(alpha_t, alpha_t-1)
// used in psi particle filter
void ung_bsm::smoother_ccov(const arma::vec& approx_y, const arma::vec& HH, 
  arma::mat& at, arma::cube& Pt, arma::cube& ccov) const {
  
  at.col(0) = a1;
  Pt.slice(0) = P1;
  arma::vec vt(n);
  arma::vec Ft(n);
  arma::mat Kt(m, n);
  
  
  for (unsigned int t = 0; t < (n - 1); t++) {
    Ft(t) = arma::as_scalar(Z.col(t * Ztv).t() * Pt.slice(t) * Z.col(t * Ztv) + HH(t));
    Kt.col(t) = Pt.slice(t) * Z.col(t * Ztv) / Ft(t);
    vt(t) = arma::as_scalar(approx_y(t) - D(t * Dtv) - Z.col(t * Ztv).t() * at.col(t));
    at.col(t + 1) = C.col(t * Ctv) + T.slice(t * Ttv) * (at.col(t) + Kt.col(t) * vt(t));
    Pt.slice(t + 1) = arma::symmatu(T.slice(t * Ttv) * (Pt.slice(t) -
      Kt.col(t) * Kt.col(t).t() * Ft(t)) * T.slice(t * Ttv).t() + RR.slice(t * Rtv));
    
    ccov.slice(t) = Pt.slice(t+1); //store for smoothing;
  }
  unsigned int t = n - 1;
  Ft(t) = arma::as_scalar(Z.col(t * Ztv).t() * Pt.slice(t) * Z.col(t * Ztv) + HH(t));
  vt(t) = arma::as_scalar(approx_y(t) - D(t * Dtv) - Z.col(t * Ztv).t() * at.col(t));
  Kt.col(t) = Pt.slice(t) * Z.col(t * Ztv) / Ft(t);
  ccov.slice(t) = arma::symmatu(T.slice(t * Ttv) * (Pt.slice(t) -
    Kt.col(t) * Kt.col(t).t() * Ft(t)) * T.slice(t * Ttv).t() + RR.slice(t * Rtv));
  
  
  arma::vec rt(m, arma::fill::zeros);
  arma::mat Nt(m, m, arma::fill::zeros);
  
  for (int t = (n - 1); t >= 0; t--) {
    arma::mat L = T.slice(t * Ttv) * (arma::eye(m, m) - Kt.col(t) * Z.col(t * Ztv).t());
    //P[t+1] stored to ccov_t
    ccov.slice(t) = Pt.slice(t) * L.t() * (arma::eye(m, m) - Nt * ccov.slice(t));
    rt = Z.col(t * Ztv) / Ft(t) * vt(t) + L.t() * rt;
    Nt = arma::symmatu(Z.col(t * Ztv) * Z.col(t * Ztv).t() / Ft(t) + L.t() * Nt * L);
    
    at.col(t) += Pt.slice(t) * rt;
    Pt.slice(t) -= arma::symmatu(Pt.slice(t) * Nt * Pt.slice(t));
  }
}

double ung_bsm::psi_filter(const arma::vec& approx_y, const arma::vec& approx_var_y,
  const arma::vec& scales, const unsigned int nsim, arma::cube& alpha, 
  arma::mat& weights, arma::umat& indices) {
  
  arma::mat alphahat(m, n);
  arma::cube Vt(m, m, n);
  arma::cube Ct(m, m, n);
  smoother_ccov(approx_y, approx_var_y, alphahat, Vt, Ct);
  conditional_cov(Vt, Ct);
  
  std::normal_distribution<> normal(0.0, 1.0);
  
  
  for (unsigned int i = 0; i < nsim; i++) {
    arma::vec um(m);
    for(unsigned int j = 0; j < m; j++) {
      um(j) = normal(engine);
    }
    alpha.slice(i).col(0) = alphahat.col(0) + Vt.slice(0) * um;
  }
  
  std::uniform_real_distribution<> unif(0.0, 1.0);
  arma::vec normalized_weights(nsim);
  double loglik = 0.0;
  
  weights.col(0) = exp(log_weights(approx_y, approx_var_y, 0, alpha) - scales(0));
//weights.col(0) = exp(weights.col(0) - weights.col(0).max());
  double sum_weights = arma::sum(weights.col(0));
  if(sum_weights > 0.0){
    normalized_weights = weights.col(0) / sum_weights;
  } else {
    return -arma::datum::inf;
  }
  loglik = log(sum_weights / nsim);
  
  
  for (unsigned int t = 0; t < (n - 1); t++) {
    arma::vec r(nsim);
    for (unsigned int i = 0; i < nsim; i++) {
      r(i) = unif(engine);
    }
    indices.col(t) = stratified_sample(normalized_weights, r, nsim);
    
    arma::mat alphatmp(m, nsim);
    
    for (unsigned int i = 0; i < nsim; i++) {
      alphatmp.col(i) = alpha.slice(indices(i, t)).col(t);
    }
    for (unsigned int i = 0; i < nsim; i++) {
      arma::vec um(m);
      for(unsigned int j = 0; j < m; j++) {
        um(j) = normal(engine);
      }
      alpha.slice(i).col(t + 1) = alphahat.col(t + 1) +
        Ct.slice(t + 1) * (alphatmp.col(i) - alphahat.col(t)) + Vt.slice(t + 1) * um;
    }
    
    weights.col(t + 1) = exp(log_weights(approx_y, approx_var_y, t + 1, alpha) - scales(t+1));
    //weights.col(t + 1) = exp(weights.col(t + 1) - weights.col(t + 1).max());
    double sum_weights = arma::sum(weights.col(t + 1));
    if(sum_weights > 0.0){
      normalized_weights = weights.col(t + 1) / sum_weights;
    } else {
      return -arma::datum::inf;
    }
    loglik += log(sum_weights / nsim);
    
  }
  return loglik;
}
