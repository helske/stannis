#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4bsm_approx_mod) {


    class_<rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> >("model_bsm_approx")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_bsm_approx_namespace::model_bsm_approx, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4ll_approx_mod) {


    class_<rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> >("model_ll_approx")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_ll_approx_namespace::model_ll_approx, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4llt_approx_mod) {


    class_<rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> >("model_llt_approx")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_llt_approx_namespace::model_llt_approx, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4poisson_ll_mod) {


    class_<rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> >("model_poisson_ll")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_poisson_ll_namespace::model_poisson_ll, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4x_bsm_approx_mod) {


    class_<rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> >("model_x_bsm_approx")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_x_bsm_approx_namespace::model_x_bsm_approx, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4x_llt_approx_mod) {


    class_<rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> >("model_x_llt_approx")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_x_llt_approx_namespace::model_x_llt_approx, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4x_llt_poisson_mod) {


    class_<rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> >("model_x_llt_poisson")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_x_llt_poisson_namespace::model_x_llt_poisson, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
