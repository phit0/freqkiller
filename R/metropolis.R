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
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
#' \dontrun{
#' # default arguments
#' data("PlantGrowth")
#' attach(PlantGrowth)
#' mh1 <- frequentistkiller(weight ~ group,
#'                         dist = "normal", number_it = 2000)
#' summary(mh1)
#' matplot(mh1$chain, type = "l", col = seq(1,4), ylim = c(-1, 6))
#' legend("center", legend = colnames(mh1$chain), col = seq(1,4),lty = 1)
#'
#' # poisson distributed data with default parameters
#'
#' set.seed(42)
#' n <- 100
#' x <- rnorm(n)
#' beta <- c(4, 1.2)
#' lambda <- exp(beta[1] + x*beta[2])
#' (y <-  rpois(n, lambda))
#'
#' mh2 <- frequentistkiller(y ~ x, dist = "poisson", number_it = 2000)
#' summary(mh2)
#' }
frequentistkiller <- function(formula, dist, beta_start = "ml_estimate",
                     a0 = 0.001, b0 = 0.001, number_it, m = beta_start,
                     M = diag(ncol(model.matrix(formula))), thinning_lag = 1, burnin = 500){
  # check if default or manual startvalue
  if(is.character(beta_start)) {
    if (beta_start == "ml_estimate") {
    beta_start <- beta_init(formula,dist)
    }else {
    stop("beta_start can be either a numeric
               vector of appropriate length or default \"ml_estimator\"")
      }
  }

  # increase number of iterations if thinning will be performed
  num_it <- (number_it + burnin) * thinning_lag

  if (dist == "poisson") {

    # run algorithm
    chain <- metroPois(formula, beta_start, m, M, num_it, dist)
  }
  else if  (dist == "normal") {
    chain <- metroNorm(formula, beta_start, a0, b0, m, M, num_it, thinning_lag, dist)
  }
  else if (dist == "bernoulli"){
    chain <- metroBer(formula, beta_start, m, M, num_it, dist)
  }
  else {
    stop("Wrong distribution name. Choose one of the implemented distributions:
         \"normal\", \"poisson\" or \"bernoulli\".")
  }
  # cut off burn in phase
  chain <- chain[burnin:num_it, ]
  # thinning if desired
  if (thinning_lag > 1) {
    chain <- chain[seq(1, nrow(chain), thinning_lag), ]
  }

  # gather objects for the output in a list
  result <- list(chain = chain, thinning_lag = thinning_lag, number_it = number_it,
                 beta_start = beta_start, m = m, M = M,
                 burnin = burnin, dist = dist, formula = formula)
  # define a second class "metrohas" for the output, in order to use the summary()
  class(result) <- append("metrohas", "list")

  print("DONE")
  return(result)
}

