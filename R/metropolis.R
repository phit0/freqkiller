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
metrohas <- function(formula, data, startvalue, anzahl_sim){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  chain <- array(dim = c(anzahl_sim + 1, length(startvalue)))

  chain[1,] <- startvalue

  for (i in 1:anzahl_sim) {
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

hey = "saludos"

hallo = "A buenahora companiero"
