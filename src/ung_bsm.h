#ifndef UNG_BSM_H
#define UNG_BSM_H

#include <RcppArmadillo.h>

class ung_bsm {

public:

  ung_bsm(const Rcpp::List& model, const unsigned int seed);
  
  // update model given the parameters theta
  void set_theta(const arma::vec& theta);
  // extract theta from the model
  arma::vec get_theta() const;
  // compute covariance matrices RR and regression part
  void compute_RR();
  void compute_xbeta() { xbeta = xreg * beta; }

  arma::vec log_obs_density(const unsigned int t, 
    const arma::cube& alpha) const;
  // bootstrap filter  
  double bsf_filter(const unsigned int nsim, arma::cube& alphasim, 
    arma::mat& weights, arma::umat& indices);
  
  arma::vec y;
  arma::mat Z;
  arma::cube T;
  arma::cube R;
  arma::cube Q;
  arma::vec a1;
  arma::mat P1;
  arma::mat xreg;
  arma::vec beta;
  arma::vec D;
  arma::mat C;
  
  const unsigned int Ztv;
  const unsigned int Ttv;
  const unsigned int Rtv;
  const unsigned int Dtv;
  const unsigned int Ctv;
  
  const unsigned int n;
  const unsigned int m;
  const unsigned int k;
  
  arma::cube RR;
  arma::vec xbeta;
  
  std::mt19937 engine;
  
  double phi;
  arma::vec u;
  unsigned int distribution;
  bool phi_est;

private:
  const bool slope;
  const bool seasonal;
  const bool noise;
  const arma::uvec fixed;
  const bool level_est;
  const bool slope_est;
  const bool seasonal_est;
  unsigned int seed;
};

#endif
