#'
#'
#'
#' @param formula object of the type formula or a one that can be coerced into the class formula
#' @param beta_start a beta vector
#' @param sigma2_start starting value for the variance, here only necessary for normal distributed data
#' @param a0 starting value for shape (a) parameter of the inverse gamma prior for \Sigma^2
#' @param b0 starting value for scale(b) parameter of the inverse gamma prior for  \Sigma^2
#' @param anzahl_sim number n of betas to sample, determines the length of the chain
#' @param m prior mean of beta
#' @param M prior covariance matrix of beta
#' @param thinning_lag if dist = normal, integer > 0 to choose lag  to thin the chain of sigmas
#'
#' @return chain of n betas that will be ploted as a histogram or line plot to show the sampling path
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
metrohas <- function(formula, dist, sigma2_start = 1, beta_start,
                     a0 = 0.001, b0 = 0.0001, anzahl_sim, m = rep(0,length(beta_start)),
                     M = diag(length(beta_start)), thinning_lag = 0){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  # check dimensions are ok
  if (dim(X)[2] != length(beta_start)) {
    stop("X and beta_start have non-conformable dimensions")
  }
  if (length(y) != dim(X)[1]){
    stop("response and covariables are of different length")
  }
  if (length(sigma2_start) != 1) {
    stop("sigma2_start should be scalar")
  }
  if (length(m) != length(beta_start)) {
    stop("m should be a vector of length ", ncol(X))
  }
  if (ncol(M) != length(beta_start)){
    stop("M should be a square matrix with dimensions ", ncol(X), "x", ncol(X))
  }

  # check nonsingularity of M
  if (nrow(M) != ncol(M) | qr(M)$rank != ncol(M)) {
    stop("M must be a full rank square matrix")
  }else if (det(M) <= 0) {
    stop("M must be positive definite")
    }
  #######################################################
                    # Choose  distribution
  #######################################################
  if (dist == "poisson") {
    result <- metroPois(formula, beta_start, m, M, anzahl_sim, dist)
  }
  else if  (dist == "normal") {
    result <- metroNorm(formula, beta_start, sigma2_start, a0, b0, m, M, anzahl_sim, thinning_lag, dist)
  }
  else if (dist == "bernoulli"){
    result <- metroBer(formula, beta_start, m, M, anzahl_sim, dist)
  }
  else {
    stop("Wrong distribution name. Choose one of the implemented distributions:
         \"normal\", \"poisson\" or \"bernoulli\".")
  }

  print("DONE")
  return(result)
}

