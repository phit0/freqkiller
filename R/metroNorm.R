#############################################
###         MCMC for normal data          ###
#############################################
metroNorm <- function(y, X, beta_start, a0, b0, m, M, number_it, dist){

  M_1 <- solve(M)
  M_det <- det(M)

  # matrix for betas
  chain <- matrix(NA, nrow = number_it + 1, ncol = length(beta_start))
  chain[1,] <- beta_start

  # vector for sigmas
  s_chain <- matrix(NA, nrow = number_it + 1, ncol = 1)

  # initializing sigma2 chain
  # a is fixed
  a_t <- length(y) / 2 + a0
  # updating first b
  b_t <- b_func(chain[1, ], y, X, b0)
  # sampling first sigma
  sigma2_t <- sigma_gibbs(a_t, b_t)
  # adding to the chain
  s_chain[1] <- sigma2_t

  for (i in 1:number_it) {

    # IWLS; beta_t is only a placeholder here
    F_t <- fisher_func(sigma2_t, beta_t, y, X, M_1, dist)
    mu_t <- mu_func(F_t, sigma2_t, beta_t, y, X, M_1, m, dist)

    # Picking proposal
    proposal <- proposalfunction(mu_t, chol2inv(chol(F_t)))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal, y, X, M_1, dist)
    mu_star <- mu_func(F_star, sigma2_t, proposal, y, X, M_1, m, dist)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)

    # Prior and log-likelihood functions
    prior_t <- prior_func(chain[i, ], m, M_1, M_det)
    prior_star <- prior_func(proposal, m, M_1, M_det)
    loglik_t <- loglik_func(chain[i, ], sigma2_t, y, X, dist)
    loglik_star <- loglik_func(proposal, sigma2_t, y, X, dist)

    # Caltulating the logarithmized acceptance probability
    alpha <- min(c(prior_star + loglik_star + q_cond_star -
                     prior_t - loglik_t - q_cond_t, 0))
    # Check if alpha is NaN
    if (any(is.nan(alpha))) {
      stop("Entries of the Fisher matrix are infinite, try with different
           starting values...")
    }

    # Sampling decision to the chain
    if (log(runif(1)) < alpha) {
      chain[i + 1, ] <- proposal
    }else{
      chain[i + 1, ] <- chain[i, ]
    }

    # updating sigma2 using a Gibbs sampler
    b_t <- b_func(chain[i, ], y, X, b0)
    sigma2_t <- sigma_gibbs(a_t, b_t)
    s_chain[i + 1] <- sigma2_t

  } # End of the loop

  # Warning message for the user if the proposals were not accepted
  if (all(chain[1:10, 1] == chain[number_it - 10:number_it, 1])) {
    warning("Proposals were apparently not accepted in the chain...
            Try different starting values or use the default \"ml_estimate\".")
  }

  # adding covariable names
  colnames(chain) <- colnames(X)
  colnames(s_chain) <- "variance"

return(cbind(chain, s_chain))
}
