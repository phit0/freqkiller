
metroNorm <- function(sigma2_start, beta_start, a0, b0, anzahl_sim){

  chain <- matrix(NA, nrow = anzahl_sim + 1, ncol = length(beta_start))
  chain[1,] <- beta_start

  s_chain <- array(dim = anzahl_sim + 1)
  s_chain[1] <- sigma2_start
  sigma2_t <- sigma2_start

  a_t <- a0
  b_t <- b0

  for (i in 1:anzahl_sim) {

    # IWLS
    F_t <- fisher_func(sigma2_t, beta_t)
    mu_t <- mu_func(sigma2_t, beta_t)

    # Pick proposal
    proposal <- proposalfunction(mu_func(sigma2_t,  beta_t), sigma = solve(fisher_func(sigma2_t,beta_t)))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal)
    mu_star <- mu_func(sigma2_t, proposal)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)  #invert in function to avoid re-inverting

    #Posterior
    prior_t <- prior_func(chain[i,])
    prior_star <- prior_func(proposal)

    loglik_t <- loglik_func(chain[i,], sigma2_t)
    loglik_star <- loglik_func(proposal, sigma2_t)

    alpha <- min(c((prior_star + loglik_star + q_cond_star)
                   / (prior_t + loglik_t + q_cond_t), 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }

    # update sigma
    a_t <- a_func(y, a_t)
    b_t <- b_func(y, chain[i + 1, ], b_t)
    sigma2_t <- sigma_gibbs(a_t, b_t)
    s_chain[i+1] <- sigma2_t
  }

return(data.frame(chain, sigma2 = s_chain))
}
