#############################################
###         MCMC for bernoulli data       ###
#############################################
metroBer <- function(y, X, beta_start, m, M, number_it, dist, notify){

  M_1 <- solve(M)
  M_det <- det(M)

  # matrix for betas
  chain <- matrix(NA, nrow = number_it + 1, ncol = length(beta_start))
  chain[1,] <- beta_start

  for (i in 1:number_it) {

    beta_t <- chain[i,]
    # This holds due to var(y) = h(eta)*(1-h(eta)) = dh(eta)
    sigma2_t <- dh(X %*% beta_t)

    # IWLS
    F_t <- fisher_func(sigma2_t, beta_t, y, X, M_1, dist)
    mu_t <- mu_func(sigma2_t, beta_t, y, X, M_1, m, dist)

    # Sampling random proposal
    proposal <- proposalfunction(mu_t, solve(F_t))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal, y, X, M_1, dist)
    mu_star <- mu_func(sigma2_t, proposal, y, X, M_1, m, dist)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)

    # Prior and log-likelihood functions
    prior_t <- prior_func(chain[i,], m, M_1, M_det)
    prior_star <- prior_func(proposal, m, M_1, M_det)
    loglik_t <- loglik_func(chain[i,], sigma2_t, y, X, dist)
    loglik_star <- loglik_func(proposal, sigma2_t, y, X, dist)

    # Calculating the logarithmized acceptance probability
    alpha <- min(c(prior_star + loglik_star + q_cond_star -
                     prior_t - loglik_t - q_cond_t, 0))
    # Check if alpha is NaN
    if (any(is.nan(alpha))) {
      stop("alpha is NaN due to unlikeliy starting values.")
    }
    # Sampling decision to the chain
    if (log(runif(1)) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
  } # End of the loop

  # Warning message for the user if the proposals were not accepted
  if (all(chain[1:10, 1] == chain[number_it - 10:number_it, 1])) {
    warning("Proposals were apparently not accepted in the chain...
                Try different starting values or use the default \"ml_estimate\".")
  }


  # adding covariable names
  colnames(chain) <- colnames(X)
  return(chain)
}
