#############################################
###         MCMC for normal data          ###
#############################################
metroNorm <- function(formula, beta_start, a0, b0, m, M, number_it, dist){

  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  M_1 <- solve(M)
  M_det <- det(M)

  # matrix for betas
  chain <- matrix(NA, nrow = number_it + 1, ncol = length(beta_start))
  chain[1,] <- beta_start
  # vector for alphas
  # alphas <- matrix(data = NA, nrow = number_it + 1)

  # vector for sigmas
  s_chain <- matrix(NA, nrow = number_it + 1, ncol = 1)

  # initialize sigma2 chain
  a_t <- length(y)/2 + a0 # a is fixed
  b_t <- b_func(chain[1, ], y, X, b0) # update first b
  sigma2_t <- sigma_gibbs(a_t, b_t) # sample first sigma
  s_chain[1] <- sigma2_t # add to chain

  for (i in 1:number_it) {

    # IWLS
    F_t <- fisher_func(sigma2_t, beta_t, y, X, M_1, dist)
    mu_t <- mu_func(sigma2_t, beta_t, y, X, M_1, m, dist)

    # Pick proposal
    proposal <- proposalfunction(mu_t, solve(F_t))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal, y, X, M_1, dist)
    mu_star <- mu_func(sigma2_t, proposal, y, X, M_1, m, dist)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)  #invert in function to avoid re-inverting

    # Posterior
    prior_t <- prior_func(chain[i, ], m, M_1, M_det)
    prior_star <- prior_func(proposal, m, M_1, M_det)

    # likelihoods
    loglik_t <- loglik_func(chain[i, ], sigma2_t, y, X, dist)
    loglik_star <- loglik_func(proposal, sigma2_t, y, X, dist)

    # acceptance probability
    alpha <- min(c(prior_star + loglik_star + q_cond_star - prior_t - loglik_t - q_cond_t, 0))
    # add alphas to output
    # alphas[i] <- alpha

    if (log(runif(1)) < alpha) {
      chain[i + 1, ] <- proposal
    }else{
      chain[i + 1, ] <- chain[i, ]
    }

    # update sigma2
    b_t <- b_func(chain[i, ], y, X, b0)
    sigma2_t <- sigma_gibbs(a_t, b_t)
    s_chain[i + 1] <- sigma2_t

    # Check for startvalue issue at iteration 1000
    if (i == 1000 & length(unique(chain[1:1000, ])) == 2) {
      warning("Proposals are not being accepted in the chain...
                Try different starting values or use the default \"ml_estimate\".")
    }
  }

  # add covariable names
  colnames(chain) <- colnames(X)
  colnames(s_chain) <- "variance"

return(cbind(chain, s_chain))
}
