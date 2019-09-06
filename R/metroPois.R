#############################################
###         MCMC for poisson data         ###
#############################################
metroPois <- function(formula, beta_start, m, M, number_it, dist){

  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  M_1 <- solve(M)
  M_det <- det(M)
  # matrix for betas
  chain <- matrix(NA, nrow = number_it + 1, ncol = length(beta_start))
  chain[1,] <- beta_start
  # vector for alphas
  alphas <- matrix(data = NA, nrow = number_it + 1)

  for (i in 1:number_it) {

    beta_t <- chain[i,]

    lambda_t <- exp(X%*%beta_t)

    # IWLS
    F_t <- fisher_func(lambda_t, beta_t, y, X, M_1, dist)
    mu_t <- mu_func(lambda_t, beta_t, y, X, M_1, m, dist)

    # Pick proposal
    proposal <- proposalfunction(mu_t, solve(F_t))

    # IWLS
    F_star <- fisher_func(lambda_t, proposal, y, X, M_1, dist)
    mu_star <- mu_func(lambda_t, proposal, y, X, M_1, m, dist)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)  #invert in function to avoid re-inverting

    #Posterior
    prior_t <- prior_func(chain[i,], m, M_1, M_det)
    prior_star <- prior_func(proposal, m, M_1, M_det)

    loglik_t <- loglik_func(chain[i,], lambda_t, y, X, dist)
    loglik_star <- loglik_func(proposal, lambda_t, y, X, dist)

    alpha <- min(c(prior_star + loglik_star + q_cond_star - prior_t - loglik_t - q_cond_t , 0))
    # add alphas to output
    alphas[i] <- alpha

    if (log(runif(1)) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }

  }
  return(data.frame(chain, alpha = alphas))
}
