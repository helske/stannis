#' Importance Sampling Type Correction of Markov Chain Monte Carlo with Stan
#'
#' The package provides an importance sampling type corrected Markov 
#' chain Monte Carlo (MCMC) algorithm for structural Poisson time series, 
#' where the MCMC is performed using Hamiltonian Monte Carlo provided by Stan, 
#' and the correction is performed with particle filtering. The 
#' package also contains a vignette providing a small simulation study which compares 
#' variants of random-walk pseudo-marginal MCMC with HMC, and the IS-corrected versions of these.
#'
#' @docType package
#' @name stannis
#' @aliases stannis
#' @importFrom Rcpp evalCpp
#' @useDynLib stannis, .registration = TRUE
NULL
