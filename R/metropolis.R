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
#' log-prior distribution and logarithmized conditional proposal density
#' compared to the previous beta. If the new value performs good enough,
#' it has a higher chance to be stored in the chain of betas. If rejected, the
#' current beta is stored in the chain and the algorithm proceeds to the next
#' iteration.
#'
#' @param formula object of the type formula or a one that can be coerced into
#' the class formula.
#' @param data \code{\link{data.frame}} or \code{\link{matrix}} with column names
#' that correspont to the values specified in the formula.
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
#' @param thinning_lag integer from \eqn{>0}
#' @param burnin integer that indicates how many samples will be cut off at the
#' beginnign of the chain, the default is 500.
#' @param notify If \code{notify = FALSE}, the user will not receive any
#' notifications when the algoritm runs.
#'
#' @details \code{dist} specified by the user as "bernoulli", "normal" or
#' "poisson" according to the assumptions about the response variable.
#' If \code{dist = "normal"}, a hybrid algorithm is executed, where samples for
#' \eqn{\sigma^2} are obtained with a Gibbs sampler updating the full
#' conditionals of a conjugate inverse gamma prior.
#' By default, the \code{beta_start} are set to the maximum likelihood estimator
#' for the regression model as estimated by \code{glm()}. \cr \cr
#' The starting values \code{a0, b0} are only needed if \code{dist = "normal"},
#' as the variance cannot be expressed by means of the expected value, as it is
#' with the other two distributions, that have only one dependent parameter.
#' Note that it might be necessary to increase the amount of iterations fixed
#' by \code{number_it}, if a larger burn in phase is selected or
#' \code{thinning_lag} is specified in order to reach a suficcient amount of
#' samples from the posterior. The remaining amount of samples are calculated
#' via \cr \code{n = (number_it - burnin) / (thinning_lag)}. \cr
#' The prior distribution of \eqn{\beta} has expectation \code{m} with
#' \code{beta_start} as default and covariance matrix \code{M}. The default is
#' the identity Matrix.\cr
#' \code{thinning_lag} specifies at what intervall the obtained chain will be
#' thinned out. If \code{thinning_lag = 1} every element of the chain will be
#' returned. If \code{thinning_lag = 10} obly every \code{10th} element will be
#' returned, and so on.
#'
#' @return A list  of class \code{\link{frequentistkiller}} with
#' the following elements:  \cr
#' \itemize{
#'  \item \code{chain}: A vector (if univariate without intercept) or dataframe
#'  with the sampled values for each element of \eqn{\beta} in the respective
#'  column.
#'  \item  The values of all arguments of the function call as specified by the
#'  user}
#'
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#' \dontrun{
#' ## Normal distribution
#' # all function arguments explicit
#' data("PlantGrowth")
#' mh1 <- frequentistkiller(weight ~ group, data = PlantGrowth, dist = "normal",
#'                          number_it = 10000, beta_start = c(3, 0, 0),
#'                          m = c(1, 0, 0), a0 = 0.001, b0 = 0.001,
#'                          M = 10*diag(c(1,1,1)), thinning_lag = 10,
#'                          burnin = 100)
#' summary(mh1)
#' matplot(mh1$chain, type = "l", col = seq(1,4), ylim = c(-1, 6),
#' main = "Traceplot")
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
#' mh2 <- frequentistkiller(y ~ x, dist = "poisson", number_it = 2000)
#' summary(mh2)
#'
#' ## Bernoulli distributed data
#' n <- 100
#' a <- runif(n)
#' df3 <- data.frame(y = rbinom(n, 1, exp(eta)/(1 + exp(eta))), x = a)
#' mh3 <- frequentistkiller(y ~ x, df3, dist = "bernoulli", number_it = 1000,
#'                                burnin = 0)
#' summary(mh3)
#' }
frequentistkiller <-
  function(formula, data = NULL, dist, beta_start = "ml_estimate",
           a0 = 0.001, b0 = 0.001, number_it, m = rep(0, length(beta_start)),
           M = diag(length(beta_start)), thinning_lag = 1, burnin = 100,
           notify = TRUE) {

  mf <- model.frame(formula, data)

  # checking if response is numeric
  if(!is.numeric(model.response(mf, type = "any"))) {
    stop("response variable has to be numeric")
  }
  y <- model.response(mf, "numeric")
  X <- model.matrix(mf, data)
  p <- ncol(X)

  # checking if default or manual startvalue
  if(is.character(beta_start)) {
    if (beta_start == "ml_estimate") {
    beta_start <- beta_init(mf, dist)
    } else {
    stop("beta_start can be either a numeric
               vector of appropriate length or default \"ml_estimate\"")
      }
  }
  # checking dimension of beta_start
  if (length(beta_start) != p) {
    stop("length of \"beta_start\" must equal the number of covariables")
  }
  # checking if positive and numeric
  if (!(is.numeric(number_it) & is.numeric(burnin) & is.numeric(thinning_lag)
        & (number_it > 0) & (burnin >= 0) & (thinning_lag > 0))) {
    stop("\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  }
  # checking if integer and scalar
  if (!((round(number_it) == number_it) & (round(burnin) == burnin)
      & (round(thinning_lag) == thinning_lag) & (length(number_it) == 1L)
      & (length(burnin) == 1L) & (length(thinning_lag) == 1L))) {
    stop("\"number_it\" \"burnin\" and \"thinning_lag\"
         have to be scalar integers")
  }
  # checking if burnin and thinning do not exceed number of iterations
  if ((number_it <= burnin) | (number_it <= thinning_lag)) {
    stop("\"number_it\" must be larger than \"burnin\" and \"thinning_lag\"")
  }
  # checking if a0 and b0 are defined properly
  if (!((is.numeric(a0)) & (is.numeric(b0))
        & (length(a0) == 1L & length(b0) == 1L))){
    stop("a0 and b0 have to be either numeric or integer scalar values")
  }
  # checking dimensions of m
  if (p != length(m)) {
    stop("\"m\" must be of the same lenght as \"beta_start\"!")
  }
  # checking dimensions of M
  if (!any(dim(M) == c(p, p))) {
    stop(paste("\"M\" must be a square matrix of ", p, "x", p))
  }

  # select a distribution and run the Metropolis-Hastings algorithm
  chain <- switch(dist,
    "poisson" = metroPois(y, X, beta_start, m, M, number_it, dist, notify),
    "normal" = metroNorm(y, X, beta_start, a0, b0, m, M, number_it, dist),
    "bernoulli" =  metroBer(y, X, beta_start, m, M, number_it, dist, notify),
    stop("Wrong distribution name. Choose one of the implemented distributions:
         \"normal\", \"poisson\" or \"bernoulli\"."))

  # cuting off burn in phase
  chain <- chain[burnin:number_it, ]

  # thinning if desired
  if (thinning_lag > 1) {
    chain <- chain[seq(1, nrow(chain), thinning_lag), ]
  }

  # gathering objects for the output in a list
  result <- list(chain = chain, thinning_lag = thinning_lag, number_it = number_it,
                 beta_start = beta_start, a0 = a0, b0 = b0, m = m, M = M,
                 burnin = burnin, dist = dist, formula = formula)

  # defining a second class "frequentistkiller" for the output, in order to use
  # the summary() function
  class(result) <- append("frequentistkiller", "list")
  if (notify) {
    print("DONE")
  }

  return(result)
}

