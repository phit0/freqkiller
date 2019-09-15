#' Sample betas from IWLS proposal densities via MCMC
#'
#' @description The framework of the model are bayesian generalized linear models
#' of the form
#' \deqn{E(y) = h(X\beta)}
#' where \eqn{h} is the response function. Further, a prior of the form
#' \eqn{\beta ~ N(m, M)} is specified. \cr \cr
#' Available distributions and associated link functions are:
#' \itemize{
#'   \item "bernoulli" - logit link
#'   \item "normal" - identity link
#'   \item "poisson" - log link}
#' The function executes the Metropolis Hashtings Algorithm, where a Normal
#' proposal density based on prior information and the current beta is generated
#' (and updated) by means of iteratively weighted least squares. A random sample
#' from the proposal density is then evaluated according to the loglikelihood,
#' log-prior distribution and logarithmized conditional proposal density compared to the
#' previous beta. If the new value performs good enough, it has a higher chance
#' to be stored in the chain of betas. If rejected, the current beta is stored
#' in the chain and the algorithm proceeds to the next iteration.
#'
#' @param formula object of the type formula or a one that can be coerced into
#' the class formula.
#' @param dist character element either "bernoulli", "normal" or "poisson".
#' @param beta_start vector of appropriate length to specify the starting values
#' for the regression parameter \eqn{\beta} in the Metropolis-Hastings algorithm.
#' @param a0 scalar number, starting parameter "shape" for the iverse gamma (IG)
#' distribution of the Gibbs sampler for the variance.
#' @param b0 scalar number, starting parameter "rate" for the IG distribution.
#' @param number_it number n of iterations for the Metropolis-Hastings algorithm.
#' @param m numeric vector or single number of the same length as \eqn{beta}.
#' @param M a numeric symetric noningular square Matrix with same size as the
#' amount of regression parameters of interest.
#' @param thinning_lag integer from \eqn{[1:100]}
#' @param burnin integer that indicates how many samples will be cut off at the
#' beginnign of the chain, the default is 500.
#'
#' @details \code{dist} specified by the user as "bernoulli", "normal" or
#' "poisson" according to the assumptions about the response variable.
#' If \code{dist = "normal"}, a hybrid algorithm is executed, where samples for
#' \code{\sigma^2} are obtained with a Gibbs sampler updating the full
#' conditionals of a conjugate inverse gamma prior.
#' By default, the \code{beta_start} are set to the maximum likelihood estimator
#' for the regression model as estimated by \code{glm()}. \cr \cr
#' The starting values \code{a0, b0} are only needed if \code{dist = "normal"},
#' as the variance cannot be expressed by means of the expected value, as it is
#' with the other two distributions, that have only one dependent parameter.
#' Note that it might be necessary to increase the amount of iterations fixed
#' by \code{number_it}, if a larger burn in phase is selected or \code{thinning_lag}
#' is specified in order to reach a suficcient amount of samples from the posterior.
#' The remaining amount of samples are calculated via \cr
#' \code{n = (number_it - burnin) / (thinning_lag)}. \cr
#' The prior distribution of \eqn{\beta} has expectation \code{m} with \code{beta_start} as default
#' and covariance matrix \code{M}.
#' The default is the identity Matrix.\cr
#' \code{thinning_lag} specifies at what intervall the obtained chain will be thinned out.
#' If \code{thinning_lag = 1} every element of the chain will be returned.
#' If \code{thinning_lag = 10} obly every \code{10th} element will be returned, and so on.
#'
#'
#'
#' @return A list  of class \code{"\link[=metrohas.object]{metrohas}"} with
#' the following elements:  \cr
#' \itemize{
#'  \item \code{chain}: A vector (if univariate without intercept) or dataframe
#'  with the sampled values for each element of \eqn{\beta} in the respective
#'  column.
#'  \item  The values of all arguments of the function call as specified by the user}
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' \dontrun{
#' ## Normal distribution
#' # all function arguments explicit
#' data("PlantGrowth")
#' attach(PlantGrowth)
#'
#' mh1 <- frequentistkiller(weight ~ group, dist = "normal", number_it = 10000,
#'                         beta_start = c(3, 0, 0), m = c(1, 0, 0), a0 = 0.001, b0 = 0.001,
#'                         M = 10*diag(c(1,1,1)), thinning_lag = 10, burnin = 500)
#' summary(mh1)
#' matplot(mh1$chain, type = "l", col = seq(1,4), ylim = c(-1, 6),  main = "Traceplot")
#' legend("center", legend = colnames(mh1$chain), col = seq(1,4),lty = 1)
#'
#' ## poisson distributed data
#' # default parameters
#' set.seed(42)
#' n <- 100
#' x <- rnorm(n)
#' beta <- c(4, 1.2)
#' lambda <- exp(beta[1] + x*beta[2])
#' (y <-  rpois(n, lambda))
#'
#' mh2 <- frequentistkiller(y ~ x, dist = "poisson", number_it = 2000)
#' summary(mh2)
#'
#' ## Bernoulli distributed data
#' }
frequentistkiller <- function(formula, dist, beta_start = "ml_estimate",
                     a0 = 0.001, b0 = 0.001, number_it, m = beta_start,
                     M = diag(ncol(model.matrix(formula))), thinning_lag = 1, burnin = 500){

  # check if default or manual startvalue
  if(is.character(beta_start)) {
    if (beta_start == "ml_estimate") {
    beta_start <- beta_init(formula, dist)
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
    chain <- metroNorm(formula, beta_start, a0, b0, m, M, num_it, dist)
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

