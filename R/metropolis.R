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
#'
#' @examples metrohas(c(1, 2, 3), 5000)
metrohas <- function(formula, beta_start, a0 = 0.001, b0 = 0.0001, anzahl_sim, m = rep(0,length(beta_start)), M = diag(length(beta_start))){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  chain <- array(dim = c(anzahl_sim + 1, length(beta_start)))

  chain[1,] <- beta_start
  a_t <- a0
  b_t <- b0

  for (i in 1:anzahl_sim) {
    eta <- X%*%chain[i,]

    # update sigma
    a_t <- a_func(y, a_t)
    b_t <- b_func(y, eta, b_t)
    sigma_t <- sigma_gibbs(y, eta, a_t, b_t)

    # IWLS
    W_t <- w_func(eta, sigma_t)

    Ft <- fisher_func()
    proposal <- proposalfunction(chain[i,], X, y)

    enumerator <- (post_beta(proposal, X, y) + cond_proposaldensity(chain[i,], proposal, X, y))
    nominator <- (post_beta(chain[i, ], X, y) + cond_proposaldensity(proposal, chain[i, ], X, y))
    alpha <- min(c(enumerator / nominator, 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  print("DONE")
  return(chain)
}

