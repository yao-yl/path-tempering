#' Adaptive Path Sampling and Continuous Tempering
#'
#' @docType package
#' @name pathtemp-package
#'
#' @importFrom stats approx
#' @importFrom graphics  abline hist lines par plot
#' @importFrom rstan  sampling extract optimizing

#' @description
#'
#'  The normalizing constant plays an important role in Bayesian computation. Adaptive path sampling
#'  iteratively reduces the gap between the proposal and the target density, and provide a reliable
#'  normalizing constant  estimation with practical diagnostic using importance sampling theory. When
#'  equip simulated tempering with a continuous temperature, path tempering enables efficient sampling
#'  from multimodal densities.
#'
#'
#'
#' @references Yuling Yao, Collin  Cademartori,Aki Vehtari,
#' Andrew Gelman Adaptive Path Sampling in Metastable Posterior Distributions

NULL
