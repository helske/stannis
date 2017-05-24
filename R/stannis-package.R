#' Bayesian Inference of State Space Models
#'
#' Bayesian inference of exponential family state space models
#' where the state dynamics are Gaussian but the observational density is either
#' Gaussian, Poisson, binomial, or negative binomial. The novel approach used by
#' the stannis package is based on the combination of fast Gaussian approximation,
#' efficient No-U-Turn sampler provided by Stan, and the parallelisable importance
#' sampling type corrected Markov chain Monte Carlo approach.
#'
#' @docType package
#' @name stannis
#' @aliases stannis
#' @importFrom Rcpp evalCpp
#' @useDynLib stannis, .registration = TRUE
NULL
