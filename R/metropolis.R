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
metrohas <- function(formula, dist, sigma2_start = 1, beta_start,
                     a0 = 0.001, b0 = 0.0001, anzahl_sim, m = rep(0,length(beta_start)),
                     M = diag(length(beta_start))){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]

  #####################################################
  ######   PROBLEM: Wenn X im global environment definiert wird, kommt es zum konflikt!!
  assign("X", X, envir = as.environment("package:BASS"))
  assign("y", y, envir = as.environment("package:BASS"))

  M_1 <- solve(M)
  M_det <- det(M)
  assign("M_1", M_1, envir = as.environment("package:BASS"))
  assign("M", M, envir = as.environment("package:BASS"))
  assign("M_det", M_det, envir = as.environment("package:BASS"))
  assign("m", m, envir = as.environment("package:BASS"))
  assign("dist", dist, envir = as.environment("package:BASS"))

  #environment(X) <- as.environment("package:BASS") # for functions
  if (dist == "poisson") {

    # run algorithm
    result <- metroPois(beta_start, anzahl_sim)

  }else if (dist == "normal") {

    result <- metroNorm(sigma2_start, beta_start, a0, b0, anzahl_sim)
  }


  print("DONE")
  return(result)
}

