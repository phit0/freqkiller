#############################################
###         MCMC for normal data          ###
#############################################
metroNorm <- function(formula, beta_start, sigma2_start, a0, b0, m, M, number_it, thinning_lag, dist){
  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]
  M_1 <- solve(M)
  M_det <- det(M)

  # matrix for betas
  chain <- matrix(NA, nrow = number_it + 1, ncol = length(beta_start))
  chain[1,] <- beta_start
  # vector for alphas
  alphas <- matrix(data = NA, nrow = number_it + 1)
  # vector for sigmas
  s_chain <- array(dim = number_it + 1)
  s_chain[1] <- sigma2_start
  sigma2_t <- sigma2_start

  a_t <- a0
  b_t <- b0

  for (i in 1:number_it) {

    # IWLS
    F_t <- fisher_func(sigma2_t, beta_t, y, X, M_1, dist)
    mu_t <- mu_func(sigma2_t, beta_t, y, X, M_1, m, dist)

    # Pick proposal
    proposal <- proposalfunction(mu_func(sigma2_t,  beta_t, y, X, M_1, m, dist),
                sigma = solve(fisher_func(sigma2_t, beta_t, y, X, M_1, dist)))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal, y, X, M_1, dist)
    mu_star <- mu_func(sigma2_t, proposal, y, X, M_1, m, dist)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)  #invert in function to avoid re-inverting

    # Posterior
    prior_t <- prior_func(chain[i,], m, M_1, M_det)
    prior_star <- prior_func(proposal, m, M_1, M_det)
    # likelihoods
    loglik_t <- loglik_func(chain[i,], sigma2_t, y, X, dist)
    loglik_star <- loglik_func(proposal, sigma2_t, y, X, dist)
    # acceptance probability
    alpha <- min(c(prior_star + loglik_star + q_cond_star - prior_t - loglik_t - q_cond_t, 0))
    # add alphas to output
    alphas[i] <- alpha

    if (log(runif(1)) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }

    # update sigma
    a_t <- a_func(y, a_t)
    b_t <- b_func(chain[i + 1, ], y, X, b_t)
    sigma2_t <- sigma_gibbs(a_t, b_t)
    s_chain[i+1] <- sigma2_t

  }
  # end of iterations
  if (thinning_lag > 0) {
    s_chain <- s_chain[seq(1, length(s_chain), thinning_lag)]
  }
  acf(s_chain, main = expression(paste("Autocorrelation of ", sigma^2)))

return(data.frame(chain, sigma2 = s_chain, alpha = alphas))
}
