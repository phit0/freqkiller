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
metrohas <- function(formula, sigma_start, beta_start, a0 = 0.001, b0 = 0.0001, anzahl_sim, m = rep(0,length(beta_start)), M = diag(length(beta_start))){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  chain <- array(dim = c(anzahl_sim + 1, length(beta_start)))
  chain[1,] <- beta_start
  eta_t <- X%*%beta_start

  s_chain <- array(dim = anzahl_sim + 1)
  s_chain[1] <- sigma_start
  sigma_t <- sigma_start


  a_t <- a0
  b_t <- b0

  for (i in 1:anzahl_sim) {


    # IWLS
    W_t <- w_func(eta_t, sigma_t)
    F_t <- fisher_func(X, W_t, M)
    y_wgl_t <- y_wgl_func(eta_t, y)
    mu_t <- mu_func(X, F_t, W_t, y_wgl_t, M, m)

    # Pick proposal
    proposal <- proposalfunction(mu_t, sigma = solve(F_t))

    #Update eta
    eta_star <- X%*%proposal

    # IWLS
    W_star <- w_func(eta_star, sigma_t)
    F_star <- fisher_func(X, W_star, M)
    y_wgl_star <- y_wgl_func(eta_star, y)
    mu_star <- mu_func(X, F_star, W_star, y_wgl_star, M, m)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, solve(F_star))
    q_cond_t <- cond_proposaldensity(proposal, mu_t, solve(F_t))

    #Posterior
    prior_t <- prior_func(chain[i,], m, M)
    prior_star <- prior_func(proposal, m, M)

    loglik_t <- loglik_func(eta_t, sigma_t, y)
    loglik_star <- loglik_func(eta_star, sigma_t, y)

    alpha <- min(c((prior_star + loglik_star + q_cond_star)
                   / (prior_t + loglik_t + q_cond_t), 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
    eta_t <- X%*%chain[i + 1,]
    # update sigma
    a_t <- a_func(y, a_t)
    b_t <- b_func(y, eta_t, b_t)
    sigma_t <- sigma_gibbs(a_t, b_t)
    s_chain[i+1] <- sigma_t

  }
  print("DONE")
  return(cbind(chain, s_chain))
}

