#' Sample betas from a proposal density using MCMC
#'
#'The function executes the Metropolis Hashtings Algorithm, where a Normal proposal density
#'based on prior information and the current beta is generated (and updated) by means of iteratively
#'weighted least squares. A random sample from the proposal density is then evaluated according to
#'the loglikelihood, log-prior distribution and conditional proposal density compared to the previous
#'beta. If the new value performs good enough, it has a higher chance to be stored in the chain of
#'betas. If rejected, the current beta is stored in the chain and the algorithm proceeds to the next iteration.
#'
#' @param startvalue a beta vector
#' @param anzahl_sim number n of betas to sample, determines the length of the chain
#' @param formula object of the type formula or a one that can be coerced into the class formula
#' @return chain of n betas that will be ploted
#'  as a histogram or line plot to show the sampling path
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom mvtnorm dmvnorm
#'
#' @examples metrohas(c(1, 2, 3), 5000)
metrohas <- function(formula, dist, sigma_start = 1, beta_start,
                     a0 = 0.001, b0 = 0.0001, anzahl_sim, m = rep(0,length(beta_start)),
                     M = diag(length(beta_start))){

  if (dist == "poisson") {

    # run algorithm
    result <- metroPois(formula, beta_start, anzahl_sim, m, M)

  }else if (dist == "normal") {

    result <- metroNorm(formula, sigma_start, beta_start, a0, b0, anzahl_sim, m, M)
  }


  print("DONE")
  return(result)
}

