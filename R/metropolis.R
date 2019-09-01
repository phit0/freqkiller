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
#' @examples
metrohas <- function(formula, dist, sigma2_start = 1, beta_start = "ml_estimator",
                     a0 = 0.001, b0 = 0.001, number_it, m = rep(0,ncol(model.matrix(formula))),
                     M = diag(ncol(model.matrix(formula))), thinning_lag = 0, burnin = 1){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]

  if(is.character(beta_start)){
  if (beta_start == "ml_estimator"){
    beta_start = beta_init(formula,dist)
  }
  }

  if (dist == "poisson") {

    # run algorithm
    chain <- metroPois(formula, beta_start, m, M, number_it, dist)
  }
  else if  (dist == "normal") {
    chain <- metroNorm(formula, beta_start, sigma2_start, a0, b0, m, M, number_it, thinning_lag, dist)
  }
  else if (dist == "bernoulli"){
    chain <- metroBer(formula, beta_start, m, M, number_it, dist)
  }
  else {
    stop("Wrong distribution name. Choose one of the implemented distributions:
         \"normal\", \"poisson\" or \"bernoulli\".")
  }
  chain <- chain[burnin:number_it, ]
  result <- list(chain = chain, thinning = thinning_lag, number_it = number_it,
                 beta_start = beta_start)
  class(result) <- append("metrohas", "list")
  print("DONE")
  return(result)
}

